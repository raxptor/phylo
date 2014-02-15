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
		BUFSIZE  = 4096 // char group * taxa limit
	};
	
	typedef unsigned char st_t;
	
	// optimization state for a tree or subtree
	
	struct cgroup_data
	{
		int count;
		st_t ostate[BUFSIZE];
		st_t pstate[BUFSIZE]; 
		st_t fstate[BUFSIZE];
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
	
	// Converts the values in the matrix into bit forms like this
	//
	// 1 = 0000001
	// 2 = 0000010
	// 3 = 0000100
	// 
	void optimize_cgroups(matrix::cgroup *source, cgroup_data *out, int maxnodes, int taxons)
	{
		out->count = source->count;
		
		if (out->count * taxons > BUFSIZE)
		{
			std::cerr << "Bump BUFSIZE in optimize.cpp to at least " << out->count * taxons << std::endl;
			exit(-1);
		}

		for (unsigned int i=0;i<source->count;i++)
		{
			for (unsigned t=0;t<taxons;t++)		
			{
				unsigned int mtxofs = i * taxons + t;
				unsigned int dataofs = i * maxnodes + t;
				out->ostate[dataofs] = mask(source->submatrix[mtxofs]);
			}
		}
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
		int fin_ancestor[1024];
		int node[1024];
		
		DPRINT("Final pass from root [" << root << "] taxons=" << taxons);
 
		for (int i=0;i<cd->count;i++)
		{
			st_t *F = &cd->fstate[maxnodes * i];
			st_t *P = &cd->pstate[maxnodes * i];
			
			fin_ancestor[0] = P[root];
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
			
			// this is a network of just 2
//			if (node[0] < taxons)
//				queue = -1;
			
			while (queue >= 0)
			{
				const int me = node[queue];
				const int fa = fin_ancestor[queue];
				const int kid0 = net[me].c1;
				const int kid1 = net[me].c2;
				
				DPRINT("Node [" << me << "] " << net[me].c0 << "," << net[me].c1 << "," << net[me].c2);
				DPRINT("Processing " << me << " (" << kid0 << "," << kid1 << ") fa=" << fa);
				
				F[me] = fa & P[me];
				if (F[me] != fa)
				{
					// change here
					if (P[kid0] & P[kid1])
					{
						const int parent_share = (P[kid0] | P[kid1]) & fa;
						F[me] = parent_share | P[me];
					}
					else
					{
						F[me] = P[me] | fa;
					}
				}
				
				queue--;
				
				if (kid0 >= taxons)
				{
					queue++;
					fin_ancestor[queue] = F[me];
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
					fin_ancestor[queue] = F[me];
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
				const int r = taxp[i];
				const int p = taxp[i+1];
			
				int f_root = P[r] & F[p];
				if (!f_root)
					f_root = P[r];
				
				DPRINT("For terminal node (" << r << ") htu=" << p << " i computed " << (int) f_root << " from " << (int)F[p] << " and " << (int)P[r]);
				
				F[r] = f_root;
			}
			
			// ANd now root
			
		}
		
		
		return 0;
	}
	
	// Fitch single character first pass
	int single_unordered_character_first_pass_calc_length(int *fpo, int maxnodes, int root, int rootHTU, cgroup_data *cd)
	{
		int sum = 0;


		// for all characters in this group
		for (int i=0;i<cd->count;i++)
		{
			const int *bp = fpo;

			// offset to the right row into the submatrix table
			st_t *prow = &cd->pstate[maxnodes * i];
			st_t *orow = &cd->ostate[maxnodes * i];

			DPRINT("Frist pass from root(" << root <<") = " << (int)prow[root] << " (o:" << (int)orow[root] << ")");

			while (bp[0] != -1)
			{
				const int n  = bp[0];
				const int c1 = bp[1];
				const int c2 = bp[2];
				
				const int a = prow[c1];
				const int b = prow[c2];
				
				DPRINT(n << " => " << c1 << ", " << c2 << "  a=" << a << " b=" << b);

				bp += 3;
				
				int v = a & b;
				if (!v)
				{
					v = a | b;
					++sum;
				}
				
				DPRINT("writing pstate[" << n << "] [char:" << i << "] = " << v << " sum=" << sum);
				prow[n] = v;
			}
			
			int rv = prow[rootHTU] & prow[root];
			if (!rv)
			{
				rv = prow[root];
				++sum;
				DPRINT("Diff at root (" << rootHTU << "/" << root << "), scor=" << sum);
			}
			
			DPRINT("offsetroot = " << (cd->count * i + root));
			DPRINT("root@" << root << " => " << rv);
			DPRINT("prow[rootHTU] = " << (int) prow[rootHTU] << " prow[root]=" << (int)prow[root]);
			prow[root] = rv;
			
		}
		
		return sum;
	}
	
	int clip_merge_dist_unordered(cgroup_data *cd, int maxnodes, int target_root, int t0, int t1)
	{
		int sum = 0;

		// for all characters in this group
		if (t1 != network::NOT_IN_NETWORK)
		{
			for (int i=0;i<cd->count;i++)
			{
				// offset to the right row into the submatrix table
				st_t *P = &cd->pstate[maxnodes * i];
				st_t *F = &cd->fstate[maxnodes * i];
				DPRINT("target tree " << t0 << "-" << t1 << " has (" << (int)F[t0] << "|" << (int)F[t1] << ")");
				DPRINT("my computed merge is " << (int)F[target_root] << "@" << target_root);
				if (!(F[target_root] & (F[t0] | F[t1])))
					++sum;
			}
		}
		else
		{
			for (int i=0;i<cd->count;i++)
			{
				// offset to the right row into the submatrix table
				st_t *O = &cd->ostate[maxnodes * i];
				st_t *F = &cd->fstate[maxnodes * i];
				if (!(F[target_root] & O[t0]))
					++sum;
			}
		}
		return sum;
	}
	
	void prepare_source_tree_root_unordered(cgroup_data *cd, int maxnodes, int s0, int s1, int new_node)
	{
		for (int i=0;i<cd->count;i++)
		{
			// offset to the right row into the submatrix table
			st_t *P = &cd->pstate[maxnodes * i];
			st_t *F = &cd->fstate[maxnodes * i];
			F[new_node] = F[s0] | F[s1];
			DPRINT("FOR source tree " << s0 << " " << s1 << " F[" << new_node << "] = " << (int)F[new_node]);
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
		for (int k=0;k<st->unordered.count;k++)
		{
			for (int j=0;;j++)
			{
				int p = st->first_pass_order[j];
				if (p == -1)
					break;
				
				// The P values are fiddled with during computation, so restore all for all characters that are
				// being used in this tree (we could be dealing with a subtree so don't touch other nodes data) 
				DPRINT("restoring pstate[" << p << "] to " << (int) st->unordered.ostate[k * st->maxnodes + p]);
				st->unordered.pstate[k * st->maxnodes + p] = st->unordered.ostate[k * st->maxnodes + p];
			}		
			
			// root isn't listed in this order list. root_htu covers case with a 2-node network where the second
			// node doesn't appear in the first_pass_order (because root has only 1 child)
			st->unordered.pstate[k * st->maxnodes + root] = st->unordered.ostate[k * st->maxnodes + root];
			st->unordered.pstate[k * st->maxnodes + root_htu] = st->unordered.ostate[k * st->maxnodes + root_htu];
			DPRINT("restoring pstate[" << root << "] to " << (int) st->unordered.ostate[k*st->maxnodes+root]);
			DPRINT("restoring pstate[" << root_htu << "] to " << (int) st->unordered.ostate[k*st->maxnodes+root_htu]);
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
	
	
	void copy(optstate *target, optstate *source)
	{
	
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
