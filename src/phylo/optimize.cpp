#include "network.h"
#include "character.h"
#include "newick.h"

#include <emmintrin.h>

#include <iostream>
#include <cstring>
#include <cstdlib>

//#define DPRINT(x) { std::cout << x << std::endl; }
#define DPRINT(x) {};

namespace optimize
{
	enum
	{
		BUFSIZE  = 65768, // char group * taxa limit
		MAX_CHARACTERS = 2048
	};
	
	#define BLOCKSIZE 16
	
	enum
	{
		MAX_PACKUNITS = MAX_CHARACTERS / BLOCKSIZE
	};
	
	typedef unsigned char st_t;
	
	// optimization state for a tree or subtree
	
	struct cgroup_data
	{
		int count;
		int taxonwidth;
		int bufsize;
		
		st_t ostate[BUFSIZE] __attribute__ ((aligned(BLOCKSIZE)));
		st_t pstate[BUFSIZE] __attribute__ ((aligned(BLOCKSIZE)));
		st_t fstate[BUFSIZE] __attribute__ ((aligned(BLOCKSIZE)));
		st_t weights[MAX_CHARACTERS] __attribute__ ((aligned(BLOCKSIZE)));
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
	void init()
	{
	
	}
	
	inline int num_units(int count)
	{
		return ((count + BLOCKSIZE - 1) / BLOCKSIZE);
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

		//		
		out->taxonwidth = num_units(out->count) * BLOCKSIZE;
		out->bufsize = out->taxonwidth * maxnodes;

		if (out->bufsize > BUFSIZE)
		{
			std::cerr << "Bump BUFSIZE in optimize.cpp to at least " << out->bufsize << std::endl;
			exit(-1);
		}
		if (source->count > MAX_CHARACTERS)
		{
			std::cerr << "Bump MAX_CHARACTERS in optimize.cpp to at least " << source->count << std::endl;
			exit(-1);
		}

		memset(out->ostate, 0x0F, out->bufsize); //out->bufsize * out->packunits);
		memset(out->pstate, 0x0F, out->bufsize); //out->bufsize * out->packunits);
		memset(out->fstate, 0x0F, out->bufsize); //out->bufsize * out->packunits);
		memset(out->weights, 0x00, out->taxonwidth);

		for (unsigned t=0;t<taxons;t++)
		{
			for (int p=0;p<out->count;p++)
			{
				unsigned int mtxofs = p * taxons + t;
				out->ostate[t * out->taxonwidth + p] = mask(source->submatrix[mtxofs]);
			}
		}
		for (int p=0;p<out->taxonwidth;p++)
			out->weights[p] = 1;
	}
	
	void copy_cgroup(cgroup_data *target, cgroup_data *source)
	{
		target->count = source->count;
		target->bufsize = source->bufsize;
		target->taxonwidth = source->taxonwidth;
		memcpy(target->ostate, source->ostate, source->bufsize);
		memcpy(target->pstate, source->pstate, source->bufsize);
		memcpy(target->fstate, source->fstate, source->bufsize);
		for (unsigned int i=0;i<source->taxonwidth;i++)
			target->weights[i] = source->weights[i];
	}
	
	optstate* create(network::data *d, bool will_copy)
	{
		optstate *st = new optstate();
		st->matrix = d->matrix;
		st->maxnodes = d->allocnodes;
		st->first_pass_order = new int[3 * d->allocnodes]; // this will be more than needed.
		st->net = new network::node[d->allocnodes];
		st->root = -1;
		
		if (!will_copy)
		{
			optimize_cgroups(&d->matrix->unordered, &st->unordered, st->maxnodes, d->mtx_taxons);
			optimize_cgroups(&d->matrix->ordered, &st->ordered, st->maxnodes, d->mtx_taxons);
		}
		return st;
	}
		
	void print_cgroup(cgroup_data *cd, int maxnodes, int taxons)
	{
		for (int t=0;t<taxons;t++)
		{
			std::cout.width(3);
			std::cout << t << " f => ";

			const int tw = cd->taxonwidth;
			
			for (int j=0;j<tw;j++)
			{
				std::cout.width(1);
				std::cout << (int)(cd->fstate[t * tw + j]);
			}
			
			std::cout << " p => ";
			for (int j=0;j<tw;j++)
			{
				std::cout.width(1);
				std::cout << (int)(cd->pstate[t * tw + j]);
			}

			std::cout << " o => ";
			for (int j=0;j<tw;j++)
			{
				std::cout.width(1);
				std::cout << (int)(cd->ostate[t * tw + j]);
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
		int ancestor[1024];
		int node[1024];
		
		DPRINT("Final pass from root [" << root << "] taxons=" << taxons);

		st_t *FF = cd->fstate;
		st_t *PP = cd->pstate;

		const int tw = cd->taxonwidth;
		
		memcpy(&FF[tw * root], &PP[tw * root], tw);

		ancestor[0] = root;
		node[0] = net[root].c1 >= 0 ? net[root].c1 : net[root].c2;
		if (node[0] < 0)
			node[0] = net[root].c0;
		
		queue = 0; 
		
		int taxp[1024];
		int taxp_count = 2; 
		taxp[0] = root;
		taxp[1] = node[0];
		
		while (queue >= 0)
		{
			const int me = node[queue];
			const int anc = ancestor[queue];
			const int kid0 = net[me].c1;
			const int kid1 = net[me].c2;

			for (int cp=0;cp<cd->taxonwidth;cp+=BLOCKSIZE)
			{
				const int blkAncestor = tw * anc + cp;
				const	int blkMe = tw * me + cp;
				const int blkKid0 = tw * kid0 + cp;
				const int blkKid1 = tw * kid1 + cp;

				__m128i ppMe = _mm_load_si128((__m128i*)&PP[blkMe]);
				__m128i fa = _mm_load_si128((__m128i*)&FF[blkAncestor]);
				__m128i pp0 = _mm_load_si128((__m128i*)&PP[blkKid0]);
				__m128i pp1 = _mm_load_si128((__m128i*)&PP[blkKid1]);
				
				__m128i fme = _mm_and_si128(ppMe, fa);

				__m128i cmpOuter = _mm_cmpeq_epi8(fme, fa);
				__m128i parent_share = _mm_and_si128(fa, _mm_or_si128(pp0, pp1));
				
				__m128i zero = _mm_setzero_si128();
				__m128i cmpInner = _mm_cmpeq_epi8(_mm_and_si128(pp0, pp1), zero);

				__m128i subopt1 = _mm_or_si128(ppMe, parent_share);
				__m128i subopt2 = _mm_or_si128(ppMe, fa);
				
				// result of inner statement
				__m128i innerResult = _mm_or_si128(_mm_andnot_si128(cmpInner, subopt1),
				                                   _mm_and_si128(cmpInner, subopt2));
	
				__m128i totalResult = _mm_or_si128(_mm_andnot_si128(cmpOuter, innerResult),
									     _mm_and_si128(cmpOuter, fme));
			
				_mm_store_si128((__m128i*)&FF[blkMe], totalResult);
			}
										
			queue--;
			
			if (kid0 >= taxons)
			{
				queue++;
				ancestor[queue] = me;
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
				ancestor[queue] = me;
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
		for (int A=0;A<taxp_count;A+=2)
		{
			const int r = tw * taxp[A];
			const int p = tw * taxp[A+1];
			for (int cp=0;cp<tw;cp+=BLOCKSIZE)
			{
				__m128i zero = _mm_setzero_si128();
				__m128i PR = _mm_load_si128((__m128i*)&PP[r + cp]);
				__m128i FP = _mm_load_si128((__m128i*)&FF[p + cp]);
				__m128i f_root = _mm_and_si128(PR, FP);
				__m128i cmp = _mm_cmpeq_epi8(f_root, zero);
				__m128i totalResult = _mm_or_si128(_mm_andnot_si128(cmp, f_root),
									     _mm_and_si128(cmp, PR));
				_mm_store_si128((__m128i*)&FF[r + cp], totalResult);
			}
		}
				
		return 0;
	}
	
	// Fitch single character first pass
	int single_unordered_character_first_pass_calc_length(int *fpo, int maxnodes, int root, int rootHTU, cgroup_data *cd)
	{
		int sum = 0;
		const int *bp = fpo;
		
		// offset to the right row into the submatrix table
		st_t *pstate = cd->pstate;

		__m128i zero = _mm_setzero_si128();
		__m128i cmp = _mm_setzero_si128();

		int n = bp[0];
		int c1 = bp[1];
		int c2 = bp[2];

		signed char tmp[MAX_CHARACTERS];
		memset(tmp, 0x00, MAX_CHARACTERS);
		
		const int tw = cd->taxonwidth;
						
		while (n != -1)
		{
			for (int cp=0;cp<cd->taxonwidth;cp+=BLOCKSIZE)
			{
				const int resaddr = n * cd->taxonwidth + cp;
				__m128i a = _mm_load_si128((__m128i*)&pstate[tw * c1 + cp]);
				__m128i b = _mm_load_si128((__m128i*)&pstate[tw * c2 + cp]);
				__m128i t = _mm_loadu_si128((__m128i*)&tmp[cp]);
				
				// interleaved here, cmp starts out at z for the first round
				// and one extra at exit
			
				__m128i v = _mm_and_si128(a, b);
				cmp = _mm_cmpeq_epi8(v, zero);
				t = _mm_add_epi8(t, cmp);
				
				__m128i alt1 = _mm_and_si128(cmp, _mm_or_si128(a, b));
				__m128i alt2 = _mm_andnot_si128(cmp, v);
				__m128i res = _mm_or_si128(alt1, alt2);
				
				_mm_store_si128((__m128i*)&pstate[resaddr], res);
				_mm_storeu_si128((__m128i*)&tmp[cp], t);
			}
			c1 = bp[3+1];
			c2 = bp[3+2];
			n = bp[3];
			bp += 3;
		}

		for (int j=0;j<cd->count;j++)
		{
			int rv = pstate[tw * rootHTU + j] & pstate[tw * root + j];
			if (!rv)
			{
				rv = pstate[tw * root + j];
				tmp[j]--;
			}
			pstate[tw * root + j] = rv;
			sum -= tmp[j] * cd->weights[j];
		}			
			
		return sum;
	}

	void ultranode(network::data *d, int node)
	{
		cgroup_data *cd = &d->opt->unordered;
		memset(&cd->fstate[cd->taxonwidth * node], 0xff, cd->taxonwidth);
	}
	
	int clip_merge_dist_unordered(cgroup_data *cd, int maxnodes, int target_root, int t0, int t1, int max)
	{
		int sum = 0;
		
		// for all characters in this group
		if (t1 != network::NOT_IN_NETWORK)
		{
			target_root *= cd->taxonwidth;
			t0 *= cd->taxonwidth;
			t1 *= cd->taxonwidth;

			__m128i zero = _mm_setzero_si128();
			__m128i grandtot = _mm_set1_epi32(0);
			const __m128i vk0 = _mm_set1_epi8(0);
			const __m128i vk1 = _mm_set1_epi16(1);   
			__m128i supertot;
			
			// in the hope for better caching
			if (t1 < t0)
				std::swap(t0, t1);

			for (int cp=0;cp<cd->taxonwidth;cp+=BLOCKSIZE)
			{
				st_t *F = &cd->fstate[cp];
				__m128i ft0 = _mm_load_si128((__m128i*)&F[t0]); 
				__m128i ft1 = _mm_load_si128((__m128i*)&F[t1]);
				__m128i ftr = _mm_load_si128((__m128i*)&F[target_root]);
				__m128i wgh = _mm_load_si128((__m128i*)(&cd->weights[cp]));
						
				__m128i rh1 = _mm_or_si128(ft0, ft1);
				__m128i tot = _mm_and_si128(ftr, rh1);
				__m128i cmp = _mm_cmpeq_epi8(tot, zero);
				__m128i subtot = _mm_and_si128(cmp, wgh);
				
				// crazy add
				__m128i vl = _mm_unpacklo_epi8(subtot, vk0);
				__m128i vh = _mm_unpackhi_epi8(subtot, vk0);
				grandtot = _mm_add_epi32(grandtot, _mm_madd_epi16(vl, vk1));
				grandtot = _mm_add_epi32(grandtot, _mm_madd_epi16(vh, vk1));
				
				supertot = _mm_add_epi32(grandtot, _mm_srli_si128(grandtot, 8));
				supertot = _mm_add_epi32(supertot, _mm_srli_si128(supertot, 4));
				sum = _mm_cvtsi128_si32(supertot);
				
				if (sum > max)
					return 100000;
			}
		}
		else
		{
			target_root *= cd->taxonwidth;
			t0 *= cd->taxonwidth;
	
			__m128i zero = _mm_setzero_si128();
			__m128i grandtot = _mm_set1_epi32(0);
			const __m128i vk0 = _mm_set1_epi8(0);       
			const __m128i vk1 = _mm_set1_epi16(1);
			__m128i supertot;
			
			for (int cp=0;cp<cd->taxonwidth;cp+=BLOCKSIZE)
			{
				st_t *O = &cd->ostate[cp];
				st_t *F = &cd->fstate[cp];
				
				__m128i wgh = _mm_load_si128((__m128i*)(&cd->weights[cp]));
				__m128i ftr = _mm_load_si128((__m128i*)&F[target_root]);
				__m128i rh1 = _mm_load_si128((__m128i*)&O[t0]);
				__m128i tot = _mm_and_si128(ftr, rh1);
				__m128i cmp = _mm_cmpeq_epi8(tot, zero);
				__m128i subtot = _mm_and_si128(cmp, wgh);

				// crazy add
				__m128i vl = _mm_unpacklo_epi8(subtot, vk0);
				__m128i vh = _mm_unpackhi_epi8(subtot, vk0);
				grandtot = _mm_add_epi32(grandtot, _mm_madd_epi16(vl, vk1));
				grandtot = _mm_add_epi32(grandtot, _mm_madd_epi16(vh, vk1));
				
				supertot = _mm_add_epi32(grandtot, _mm_srli_si128(grandtot, 8));
				supertot = _mm_add_epi32(supertot, _mm_srli_si128(supertot, 4));
				sum = _mm_cvtsi128_si32(supertot);
				
				if (sum > max)
					return 100000;
			}			
		}
		
		return sum;
	}
	
	void prepare_source_tree_root_unordered(cgroup_data *cd, int maxnodes, int s0, int s1, int new_node)
	{
		const int tw = cd->taxonwidth;
		new_node *= tw;
		s0 *= tw;
		s1 *= tw;
		
		for (int cp=0;cp<cd->taxonwidth;cp+=BLOCKSIZE)
		{
			// offset to the right row into the submatrix table
			st_t *F = &cd->fstate[cp];
			__m128i a = _mm_load_si128((__m128i*)&F[s0]);
			__m128i b = _mm_load_si128((__m128i*)&F[s1]);
			__m128i c = _mm_or_si128(a, b);
			_mm_store_si128((__m128i*)&F[new_node], c);
		}
	}
	
	void prepare_source_tree_root(network::data *d, int s0, int s1, int new_node)
	{
		prepare_source_tree_root_unordered(&d->opt->unordered, d->allocnodes, s0, s1, new_node);
	}
	
	// source tree must have have been prepared on target_root with the edge to join
	character::distance_t clip_merge_dist(network::data *d, int target_root, int t0, int t1, int maxdist)
	{
		return clip_merge_dist_unordered(&d->opt->unordered, d->allocnodes, target_root, t0, t1, maxdist);
	}

	void reset_taxons_states(cgroup_data *cd, int *first_pass_order, int taxons, int root, int root_htu)
	{
		for (int j=0;;j++)
		{
			int p = first_pass_order[j];
			if (p == -1)
				break;
			
			// The P values are fiddled with during computation, so restore all for all characters that are
			// being used in this tree (we could be dealing with a subtree so don't touch other nodes data) 
			const int ofs = p * cd->taxonwidth;
		
			if (p < taxons)
				memcpy(&cd->pstate[ofs], &cd->ostate[ofs], cd->taxonwidth);
		}		
		
		// root isn't listed in this order list. root_htu covers case with a 2-node network where the second
		// node doesn't appear in the first_pass_order (because root has only 1 child)
		const int tw = cd->taxonwidth;
		const int ofs0 = tw * root;
		const int ofs1 = tw * root_htu;
		memcpy(&cd->pstate[ofs0], &cd->ostate[ofs0], tw);
		memcpy(&cd->pstate[ofs1], &cd->ostate[ofs1], tw);
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


		// Restore taxon nodes' pstate to o-state to properly deal with multistates.
		//
		reset_taxons_states(&st->unordered, st->first_pass_order, d->mtx_taxons, root, root_htu);

		DPRINT("Root=" << root << " rootHTU=" << root_htu);
		
		const int sum0 = 0;// single_unordered_character_first_pass_calc_length(st->first_pass_order, st->maxnodes, root, root_htu, &st->ordered);
		const int sum1 = single_unordered_character_first_pass_calc_length(st->first_pass_order, st->maxnodes, root, root_htu, &st->unordered);
		
		const int sum = sum0 + sum1;
		
		DPRINT("Optimized sum = " << sum0 << "+" << sum1 << "=" << sum);
		
		if (write_final)
		{
			single_unordered_character_final_pass(root, st->maxnodes, d->mtx_taxons, st->net, &st->unordered);
		}

		return sum;
	}

	void set_weight(optstate *st, int pos, int weight)
	{
		if (pos < st->unordered.count)
		{
			st->unordered.weights[pos] = weight;
		}
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

	void copy(optstate *target, optstate *source)
	{
		copy_cgroup(&target->ordered, &source->ordered);
		copy_cgroup(&target->unordered, &source->unordered);
	}

}
