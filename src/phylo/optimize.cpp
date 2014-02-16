#include "network.h"
#include "character.h"
#include "newick.h"

#include <iostream>
#include <cstring>
#include <cstdlib>

//#define DPRINT(x) { std::cout << x << std::endl; }
#define DPRINT(x) {};

namespace optimize
{
	enum {
		BUFSIZE  = 80000, // char group * taxa limit
		MAX_CHARACTERS = 2048
	};
	
	enum {
		BLOCKSIZE = 32 // do not change
	};
	
	typedef unsigned char st_t;
	
	// optimization state for a tree or subtree
	
	struct cgroup_data
	{
		int count;
		int packunits;
		int memwidth;
		st_t ostate[BUFSIZE];
		st_t pstate[BUFSIZE]; 
		st_t fstate[BUFSIZE];
		int weights[MAX_CHARACTERS];
	};
	
	struct optstate
	{
		int root, maxnodes;
		int *first_pass_order;
		
		cgroup_data unordered;
		cgroup_data ordered;
		
		network::node *net;
		matrix::data *matrix;
	};

	inline char mask(int val)
	{
		if (val == character::UNKNOWN_CHAR_VALUE)
			return 0xff;
			
		return 1 << val;
	}

	// ----------------------- 
	// 
	//	
	int is_single[256];
	int mindist[256];

	void init()
	{
		for (int i=0;i<256;i++)
			is_single[i] = -1;
		for (int i=0;i<8;i++)
			is_single[1 << i] = i;
			
		for (int i=0;i<256;i++)
		{
			mindist[i] = 8;
			for (int a=0;a<8;a++)
			{
				for (int b=a+1;b<8;b++)
				{
					if (((1 << a) & i) && ((1 << b) & i))
					{
						int d = a - b;
						if (d < 0) d = -d;
						if (mindist[i] > d)
							mindist[i] = d;
					}
				}
			}
		}
	}
	
	inline int num_units(int count)
	{
		return ((count + BLOCKSIZE - 1) / BLOCKSIZE);
	}

	inline int mem_width(int count)
	{
		return num_units(count) * BLOCKSIZE;
	}
	
	inline int char_ofs(cgroup_data *source, int blockindex, int charofs, int taxon)
	{
		return blockindex * BLOCKSIZE + taxon * source->memwidth + charofs;
	}
	
	// Converts the values in the matrix into bit forms like this
	//
	// 1 = 0000001
	// 2 = 0000010
	// 3 = 0000100
	// 
	void optimize_cgroups(matrix::cgroup *source, cgroup_data *out, int maxnodes, int taxons)
	{
		out->count = source->count;
		
		
		// one row is
		//              T0   T1   T2   T3   T4
		// CHARBLOCK0   ABCD ABCD ABCD ABCD ABCD
		// CHARBLOCK1   EFGH EFGH EFGH EFGH EFGH
		// CHARBLOCK3   IJ00 IJ00 IJ00 IJ00 IJ00
		//
		// memwidth=    ------------------------   = BLOCKSIZE * taxons
		// memwidth * packunits = total matrix size
		
		out->memwidth = BLOCKSIZE * maxnodes;
		out->packunits = num_units(out->count);

		if (out->memwidth * out->packunits > BUFSIZE)
		{
			std::cerr << "Bump BUFSIZE in optimize.cpp to at least " << out->memwidth * maxnodes << std::endl;
			exit(-1);
		}
		if (out->packunits * BLOCKSIZE > MAX_CHARACTERS)
		{
			std::cerr << "Bump MAX_CHARACTERS in optimize.cpp to at least " << source->count << std::endl;
			exit(-1);
		}

		memset(out->ostate, 0x01, out->memwidth * out->packunits);
		memset(out->pstate, 0x01, out->memwidth * out->packunits);
		memset(out->fstate, 0x01, out->memwidth * out->packunits);
		memset(out->weights, 0x00, out->packunits * BLOCKSIZE);

		for (unsigned int i=0;i<out->packunits;i++)
		{
			int count = std::min<int>(BLOCKSIZE, out->count - i*BLOCKSIZE);
			for (unsigned t=0;t<taxons;t++)
			{
				for (int p=0;p<count;p++)
				{
					unsigned int mtxofs = (i * BLOCKSIZE + p) * taxons + t;
					out->ostate[i * out->memwidth + BLOCKSIZE * t + p] = mask(source->submatrix[mtxofs]);
				}
			}
			for (int p=0;p<BLOCKSIZE;p++)
				out->weights[i*BLOCKSIZE+p] = 1;
		}
		/*

		for (unsigned int t=0;t<taxons;t++)
		{
			for (unsigned int i=0;i<out->count;i++)
			{
				std::cout << (int)source->submatrix[i * taxons + t] << " ";
			}
			std::cout << std::endl;				
		}
		std::cout << std::endl;		
		for (unsigned int t=0;t<taxons;t++)
		{
			for (unsigned int i=0;i<out->packunits;i++)
			{
				for (unsigned int p=0;p<BLOCKSIZE;p++)
					std::cout << (int)is_single[out->ostate[i * out->memwidth + BLOCKSIZE * t + p]] << " ";

			}
			std::cout << std::endl;				
		}
		*/
	}
	
	void copy_cgroup_weights(cgroup_data *target, cgroup_data *source)
	{
		for (unsigned int i=0;i<source->count;i++)
			target->weights[i] = source->weights[i];
	}
	
	optstate* create(network::data *d)
	{
		optstate *st = new optstate();
		st->matrix = d->matrix;
		st->maxnodes = d->allocnodes;
		st->first_pass_order = new int[3 * d->allocnodes]; // this will be more than needed.
		st->net = new network::node[d->allocnodes];
		st->root = -1;
		
		optimize_cgroups(&d->matrix->unordered, &st->unordered, st->maxnodes, d->mtx_taxons);
		optimize_cgroups(&d->matrix->ordered, &st->ordered, st->maxnodes, d->mtx_taxons);
		return st;
	}
	
	char to_char(st_t x)
	{
		if (x == 0xff)
			return '?';
		else 
		{
			int v = is_single[x];
			if (v != -1)
				return '0' + v;
			else
				return 'X';
		}
	}
	
	void print_cgroup(cgroup_data *cd, int maxnodes, int taxons)
	{
		for (int t=0;t<taxons;t++)
		{
			std::cout.width(3);
			std::cout << t << " f => ";
			for (unsigned int i=0;i<cd->count;i++)
			{
				std::cout.width(3);
				std::cout << (int)(cd->fstate[i * maxnodes + t]);
			}
			std::cout << " p => ";
			for (unsigned int i=0;i<cd->count;i++)
			{
				std::cout.width(3);
				std::cout << (int)(cd->pstate[i * maxnodes + t]);
			}
			std::cout << std::endl;
		}
	}

	void print_state(optstate *st, int maxnodes, int taxons)
	{
		std::cout << std::endl;

		if (st->ordered.count > 0)
		{
			std::cout << "=== Ordered characters (" << st->unordered.count << ")" << " ===" << std::endl;
			std::cout << std::endl;
			print_cgroup(&st->ordered, maxnodes, taxons);
			std::cout << std::endl;
		}
		
		if (st->unordered.count > 0)
		{
			std::cout << "=== Unordered characters (" << st->unordered.count << ")" << " ===" << std::endl;
			std::cout << std::endl;
			print_cgroup(&st->unordered, maxnodes, taxons);
		}
	}
	
	
	// Fitch single character final pass
	int single_unordered_character_final_pass(int root, int maxnodes, int taxons, network::node *net, cgroup_data *cd)
	{
		int queue;
		int fin_ancestor[256][8];
		int node[1024];
		
		DPRINT("Final pass from root [" << root << "] taxons=" << taxons);

		for (int i=0;i<cd->packunits;i++)
		{
			st_t *F = &cd->fstate[cd->memwidth * i];
			st_t *P = &cd->pstate[cd->memwidth * i];
			
			for (int j=0;j<BLOCKSIZE;j++)
				fin_ancestor[0][j] = P[BLOCKSIZE * root + j];
				
			node[0] = net[root].c1 >= 0 ? net[root].c1 : net[root].c2;
			if (node[0] < 0)
				node[0] = net[root].c0;
			
			queue = 0; 
			
			DPRINT(" character(" << i << "), root=" << root << " fin_ancestor=" << (int)P[root]);
			DPRINT(" first node=" << node[0]);
			
			int root_htu = node[0];
			
			int taxp[1024];
			int taxp_count = 2; 
			taxp[0] = root;
			taxp[1] = node[0];
			
			while (queue >= 0)
			{
				const int me = node[queue];
				const int kid0 = net[me].c1;
				const int kid1 = net[me].c2;
				
				const int blkMe = BLOCKSIZE * me;
				const int blkKid0 = BLOCKSIZE * kid0;
				const int blkKid1 = BLOCKSIZE * kid1;
				
				for (int j=0;j<BLOCKSIZE;j++)
				{
					const int fa = fin_ancestor[queue][j];
		
					F[blkMe + j] = fa & P[blkMe + j];
					if (F[blkMe + j] != fa)
					{
						// change here
						if (P[blkKid0 + j] & P[blkKid1 + j])
						{
							const int parent_share = (P[blkKid0 + j] | P[blkKid1 + j]) & fa;
							F[blkMe + j] = parent_share | P[blkMe + j];
						}
						else
						{
							F[blkMe + j] = P[blkMe + j] | fa;
						}
					}
				}
								
				queue--;
				
				if (kid0 >= taxons)
				{
					queue++;
					for (int j=0;j<BLOCKSIZE;j++)
						fin_ancestor[queue][j] = F[blkMe + j];
					node[queue] = kid0;
				}
				else if (kid0 >= 0)
				{
					taxp[taxp_count] = kid0;				
					taxp[taxp_count+1] = me;
					taxp_count += 2;
				}
				if (kid1 >= taxons)
				{
					queue++;
					for (int j=0;j<BLOCKSIZE;j++)
						fin_ancestor[queue][j] = F[blkMe + j];
					node[queue] = kid1;
				}
				else if (kid1 >= 0)
				{
					taxp[taxp_count] = kid1;
					taxp[taxp_count+1] = me;
					taxp_count += 2;
				}
			}
			
			// Need to treat terminal nodes specially when they can have ? in them
			for (int i=0;i<taxp_count;i+=2)
			{
				const int r = BLOCKSIZE * taxp[i];
				const int p = BLOCKSIZE * taxp[i+1];
			
				for (int j=0;j<BLOCKSIZE;j++)
				{
					int f_root = P[r + j] & F[p + j];
					if (!f_root)
						f_root = P[r + j];
					F[r + j] = f_root;
				}
			}
		}
					
		return 0;
	}
	
	// Fitch single character first pass
	int single_unordered_character_first_pass_calc_length(int *fpo, int maxnodes, int root, int rootHTU, cgroup_data *cd)
	{
		int sum = 0;

		for (int i=0;i<cd->packunits;i++)
		{
			const int *bp = fpo;
			
			// offset to the right row into the submatrix table
			st_t *prow = &cd->pstate[cd->memwidth * i];
			st_t *orow = &cd->ostate[cd->memwidth * i];
			
			int subsum[BLOCKSIZE] = {0};

			while (bp[0] != -1)
			{
				const int n  = bp[0];
				const int c1 = bp[1];
				const int c2 = bp[2];
				
				for (int j=0;j<BLOCKSIZE;j++)
				{
					const int a = prow[c1*BLOCKSIZE+j];
					const int b = prow[c2*BLOCKSIZE+j];
				
					DPRINT("OFFSET is " << cd->memwidth * i << " and " << j);
					DPRINT(n << " => " << c1 << ", " << c2 << "  a=" << a << " b=" << b);
					
					int v = a & b;
					if (!v)
					{
						v = a | b;
						++subsum[j];
					}
					
					DPRINT("writing pstate[" << n << "] [char:" << i << "] = " << v << " sum=" << sum);
					prow[n*BLOCKSIZE+j] = v;
				}
				
				bp += 3;
			}
			
			for (int j=0;j<BLOCKSIZE;j++)
			{
				int rv = prow[rootHTU*BLOCKSIZE+j] & prow[root*BLOCKSIZE+j];
				if (!rv)
				{
					rv = prow[root*BLOCKSIZE+j];
					++subsum[j];
					DPRINT("Diff at root (" << rootHTU << "/" << root << "), scor=" << sum);
				}
				prow[root*BLOCKSIZE+j] = rv;
				sum += subsum[j] * cd->weights[i * BLOCKSIZE + j];
			}
		}
		
		return sum;
	}
	
	int clip_merge_dist_unordered(cgroup_data *cd, int maxnodes, int target_root, int t0, int t1)
	{
		int sum = 0;
		
		// for all characters in this group
		if (t1 != network::NOT_IN_NETWORK)
		{
			target_root *= BLOCKSIZE;
			t0 *= BLOCKSIZE;
			t1 *= BLOCKSIZE;
		
			for (int i=0;i<cd->packunits;i++)
			{
				st_t *P = &cd->pstate[cd->memwidth * i];
				st_t *F = &cd->fstate[cd->memwidth * i];
				for (int j=0;j<BLOCKSIZE;j++)
				{
					// offset to the right row into the submatrix table
					if (!(F[target_root+j] & (F[t0+j] | F[t1+j])))
					{
						sum += cd->weights[i * BLOCKSIZE + j];
					}
				}
			}
		}
		else
		{
			target_root *= BLOCKSIZE;
			t0 *= BLOCKSIZE;
			t1 *= BLOCKSIZE;
		
			for (int i=0;i<cd->packunits;i++)
			{
				st_t *O = &cd->ostate[cd->memwidth * i];
				st_t *F = &cd->fstate[cd->memwidth * i];
				for (int j=0;j<BLOCKSIZE;j++)
				{
					if (!(F[target_root+j] & O[t0+j]))
					{
						sum += cd->weights[BLOCKSIZE * i + j];
					}
				}
			}
		}
		
		return sum;
	}
	
	void prepare_source_tree_root_unordered(cgroup_data *cd, int maxnodes, int s0, int s1, int new_node)
	{
		new_node *= BLOCKSIZE;
		s0 *= BLOCKSIZE;
		s1 *= BLOCKSIZE;
		
		for (int i=0;i<cd->packunits;i++)
		{
			// offset to the right row into the submatrix table
			st_t *P = &cd->pstate[cd->memwidth * i];
			st_t *F = &cd->fstate[cd->memwidth * i];
			for (int j=0;j<BLOCKSIZE;j++)
				F[new_node+j] = F[s0+j] | F[s1+j];
		}
	}
	
	void prepare_source_tree_root(network::data *d, int s0, int s1, int new_node)
	{
		prepare_source_tree_root_unordered(&d->opt->unordered, d->allocnodes, s0, s1, new_node);
	}
	
	// source tree must have have been prepared on target_root with the edge to join
	character::distance_t clip_merge_dist(network::data *d, int target_root, int t0, int t1)
	{
		return clip_merge_dist_unordered(&d->opt->unordered, d->allocnodes, target_root, t0, t1);
	}

	int optimize_for_tree(optstate *st, network::data *d, int root, bool write_final = false)
	{
		st->root = root;
		
		int root_htu;
		
		if (write_final)
		{
			network::treeify(d, root, st->net, st->first_pass_order);
			if (st->net[root].c0 >= 0)
				root_htu = st->net[root].c0;
			else if (st->net[root].c1 >= 0)
				root_htu = st->net[root].c1;
			else
				root_htu = st->net[root].c2;
		}
		else
		{
			network::make_traverse_order(d, root, st->first_pass_order, &root_htu);
		}

		// restore to starting value (this is for the root only maybe?)
		for (int k=0;k<st->unordered.packunits;k++)
		{
			for (int j=0;;j++)
			{
				int p = st->first_pass_order[j];
				if (p == -1)
					break;
				
				// The P values are fiddled with during computation, so restore all for all characters that are
				// being used in this tree (we could be dealing with a subtree so don't touch other nodes data) 
				for (int a=0;a<BLOCKSIZE;a++)
				{
					const int ofs = k * st->unordered.memwidth + p * BLOCKSIZE + a;
					st->unordered.pstate[ofs] = st->unordered.ostate[ofs];
				}
			}		
			
			// root isn't listed in this order list. root_htu covers case with a 2-node network where the second
			// node doesn't appear in the first_pass_order (because root has only 1 child)
			for (int a=0;a<BLOCKSIZE;a++)	
			{
				const int ofs0 = k * st->unordered.memwidth + root * BLOCKSIZE + a;
				const int ofs1 = k * st->unordered.memwidth + root_htu * BLOCKSIZE + a;
				st->unordered.pstate[ofs0] = st->unordered.ostate[ofs0];
				st->unordered.pstate[ofs1] = st->unordered.ostate[ofs1];
			}
		}

		DPRINT("Root=" << root << " rootHTU=" << root_htu);
		
		const int sum0 = 0;// single_unordered_character_first_pass_calc_length(st->first_pass_order, st->maxnodes, root, root_htu, &st->ordered);
		const int sum1 = single_unordered_character_first_pass_calc_length(st->first_pass_order, st->maxnodes, root, root_htu, &st->unordered);
		const int sum = sum0 + sum1;

		DPRINT("Optimized sum = " << sum0 << "+" << sum1 << "=" << sum);
		
		if (write_final)
		{
			DPRINT("Writing final values");
			single_unordered_character_final_pass(root, st->maxnodes, d->mtx_taxons, st->net, &st->unordered);
		}

		// -- lala lala --
		return sum;
	}
	
	void set_weight(optstate *st, int pos, int weight)
	{
		if (weight < st->unordered.count)
		{
			st->unordered.weights[pos] = weight;
		}
		else
		{
			std::cerr << "Weight out of range, pos=" << pos << " < " << st->unordered.count << std::endl;
		}
	}
	
	void copy(optstate *target, optstate *source)
	{
		copy_cgroup_weights(&target->ordered, &source->ordered);
		copy_cgroup_weights(&target->unordered, &source->unordered);
	}

	void free(optstate *s)
	{
		delete [] s->net;
		delete s;
	}
		
	character::distance_t optimize(network::data *data, int root, bool write_final)
	{
		DPRINT("[optimize] - First pass run to calculate length");
		return optimize_for_tree(data->opt, data, root, write_final);
	}
}
