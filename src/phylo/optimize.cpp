#include "network.h"
#include "character.h"
#include "newick.h"

#include <emmintrin.h>

#include <iostream>
#include <cstring>
#include <cstdlib>

//#define DPRINT(x) { std::cout << x << std::endl; }
#define DPRINT(x) {};

#define USE_SIMD

namespace optimize
{
	enum {
		BUFSIZE  = 32768, // char group * taxa limit
		MAX_CHARACTERS = 1024
	};
	
	#define BLOCKSIZE 16
	
	enum {
		MAX_PACKUNITS = MAX_CHARACTERS / BLOCKSIZE
	};
	
	typedef unsigned char st_t;
	
	// optimization state for a tree or subtree
	
	struct cgroup_data
	{
		int count;
		int packunits;
		int memwidth;
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

	struct tmpsort
	{
		int val;
		int idx;
		bool operator<(const tmpsort &b) const
		{
			return val < b.val;
		} 
	};	
	
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
			std::cout << "packunits:" << out->packunits <<  " memwidth=" << out->memwidth << std::endl;
			std::cerr << "Bump BUFSIZE in optimize.cpp to at least " << out->memwidth * out->packunits << std::endl;
			exit(-1);
		}
		if (out->packunits * BLOCKSIZE > MAX_CHARACTERS)
		{
			std::cerr << "Bump MAX_CHARACTERS in optimize.cpp to at least " << source->count << std::endl;
			exit(-1);
		}
		if (out->packunits > MAX_PACKUNITS)
		{
			std::cerr << "Check calculation of MAX_PACKUNITS in optimize.cpp" << std::endl;
			exit(-1);
		}

		memset(out->ostate, 0x03, out->memwidth * out->packunits);
		memset(out->pstate, 0x03, out->memwidth * out->packunits);
		memset(out->fstate, 0x03, out->memwidth * out->packunits);
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
			for (int p=0;p<out->count;p++)
				out->weights[i*BLOCKSIZE+p] = 1;
		}
	}
	
	void copy_cgroup(cgroup_data *target, cgroup_data *source)
	{
		target->memwidth = source->memwidth;
		target->packunits = source->packunits;
		memcpy(target->ostate, source->ostate, source->memwidth * source->packunits);
		memcpy(target->pstate, source->pstate, source->memwidth * source->packunits);
		memcpy(target->fstate, source->fstate, source->memwidth * source->packunits);
		for (unsigned int i=0;i<source->packunits*BLOCKSIZE;i++)
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
			for (int i=0;i<cd->packunits;i++)
			{
				for (unsigned int j=0;j<BLOCKSIZE;j++)
				{
					std::cout.width(1);
					std::cout << (int)(cd->fstate[i * cd->memwidth + BLOCKSIZE * t + j]);
				}
			}
			std::cout << " p => ";
			for (int i=0;i<cd->packunits;i++)
			{
				for (unsigned int j=0;j<BLOCKSIZE;j++)
				{
					std::cout.width(1);
					std::cout << (int)(cd->pstate[i * cd->memwidth + BLOCKSIZE * t + j]);
				}
			}
			std::cout << " o => ";
			for (int i=0;i<cd->packunits;i++)
			{
				for (unsigned int j=0;j<BLOCKSIZE;j++)
				{
					std::cout.width(1);
					std::cout << (int)(cd->ostate[i * cd->memwidth + BLOCKSIZE * t + j]);
				}
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
	int single_unordered_character_final_pass(int root, int maxnodes, int taxons, network::node *net, cgroup_data *cd, int i)
	{
		int queue;
		int ancestor[1024];
		int node[1024];
		
		DPRINT("Final pass from root [" << root << "] taxons=" << taxons);

		st_t *F = &cd->fstate[cd->memwidth * i];
		st_t *P = &cd->pstate[cd->memwidth * i];
		
		for (int j=0;j<BLOCKSIZE;j++)
			F[BLOCKSIZE * root + j] = P[BLOCKSIZE * root + j];
			
		ancestor[0] = root;
		node[0] = net[root].c1 >= 0 ? net[root].c1 : net[root].c2;
		if (node[0] < 0)
			node[0] = net[root].c0;
		
		queue = 0; 
		
		DPRINT(" character(" << i << "), root=" << root << " fin_ancestor=" << (int)P[root]);
		DPRINT(" first node=" << node[0]);
		
		int taxp[1024];
		int taxp_count = 2; 
		taxp[0] = root;
		taxp[1] = node[0];
		
		const int in_ofs = cd->memwidth * i;
		
		while (queue >= 0)
		{
			const int me = node[queue];
			const int anc = ancestor[queue];
			
			st_t *FF = &cd->fstate[in_ofs];
			st_t *PP = &cd->pstate[in_ofs];
		
			const int blkAncestor = BLOCKSIZE * anc;
			const int blkMe = BLOCKSIZE * me;
			
			const int kid0 = net[me].c1;
			const int kid1 = net[me].c2;
			
			const int blkKid0 = BLOCKSIZE * kid0;
			const int blkKid1 = BLOCKSIZE * kid1;
			
			#if defined(USE_SIMD)
						
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


			#else
				for (int j=0;j<BLOCKSIZE;j++)
				{
					const int fa = FF[blkAncestor];
					const int pp0 = PP[blkKid0];
					const int pp1 = PP[blkKid1];
					const int ppMe = PP[blkMe];
					const int parent_share = (pp0 | pp1) & fa;
					
					const int subopt1 = ppMe | parent_share;
					const int subopt2 = ppMe | fa;
					int fme = ppMe & fa;
					
					if (fme != fa)
					{
						if (pp0 & pp1)
						{
							fme = subopt1;
						}
						else
						{
							fme = subopt2;
						}
					}
					
					FF[blkMe] = fme;
					
					FF++;
					PP++;
				}
			#endif
							
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
			const int r = BLOCKSIZE * taxp[A];
			const int p = BLOCKSIZE * taxp[A+1];

#if defined(USE_SIMD)
			__m128i zero = _mm_setzero_si128();
	
			__m128i PR = _mm_load_si128((__m128i*)&P[r]);
			__m128i FP = _mm_load_si128((__m128i*)&F[p]);
			__m128i f_root = _mm_and_si128(PR, FP);
			__m128i cmp = _mm_cmpeq_epi8(f_root, zero);
			__m128i totalResult = _mm_or_si128(_mm_andnot_si128(cmp, f_root),
								     _mm_and_si128(cmp, PR));
			_mm_store_si128((__m128i*)&F[r], totalResult);

#else
			for (int j=0;j<BLOCKSIZE;j++)
			{
				int f_root = P[r + j] & F[p + j];
				if (!f_root)
					f_root = P[r + j];
					
				// these give the multi-state characters their final values
				// at least for now this is how those ? are treated, they get
				// final values here
				F[r + j] = f_root;
			}
#endif

		}
				
		return 0;
	}

	// Fitch single character final pass
	int final_state_reoptimization(int root, int taxons, network::node *net, cgroup_data *cd, cgroup_data *ref, int i)
	{
		int queue;
		int ancestor[1024];
		int node[1024];
		
		st_t *F = &cd->fstate[cd->memwidth * i];
		st_t *P = &cd->pstate[cd->memwidth * i];
		st_t *FREF = &ref->fstate[cd->memwidth * i];
		
		for (int j=0;j<BLOCKSIZE;j++)
			F[BLOCKSIZE * root + j] = P[BLOCKSIZE * root + j];
			
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
		
			const int blkAncestor = BLOCKSIZE * anc;
			const int blkMe = BLOCKSIZE * me;
			
			const int kid0 = net[me].c1;
			const int kid1 = net[me].c2;
			
			const int blkKid0 = BLOCKSIZE * kid0;
			const int blkKid1 = BLOCKSIZE * kid1;
			
			#if defined(USE_SIMD)
			
				__m128i ppMe = _mm_load_si128((__m128i*)&P[blkMe]);
				__m128i fa = _mm_load_si128((__m128i*)&F[blkAncestor]);
				__m128i pp0 = _mm_load_si128((__m128i*)&P[blkKid0]);
				__m128i pp1 = _mm_load_si128((__m128i*)&P[blkKid1]);
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
			
				_mm_store_si128((__m128i*)&F[blkMe], totalResult);

				__m128i fref = _mm_load_si128((__m128i*)&FREF[blkMe]);
				__m128i vcmp = (__m128i)_mm_cmpeq_ps((__m128)fref, (__m128)totalResult);
				
				if ( _mm_movemask_epi8(vcmp) == 0xffff)
				{
					--queue;
					continue;
				}
			#endif

			// Add kids			
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
			const int r = BLOCKSIZE * taxp[A];
			const int p = BLOCKSIZE * taxp[A+1];

#if defined(USE_SIMD)
			__m128i zero = _mm_setzero_si128();
	
			__m128i PR = _mm_load_si128((__m128i*)&P[r]);
			__m128i FP = _mm_load_si128((__m128i*)&F[p]);
			__m128i f_root = _mm_and_si128(PR, FP);
			__m128i cmp = _mm_cmpeq_epi8(f_root, zero);
			__m128i totalResult = _mm_or_si128(_mm_andnot_si128(cmp, f_root),
								     _mm_and_si128(cmp, PR));
			_mm_store_si128((__m128i*)&F[r], totalResult);

#else
			for (int j=0;j<BLOCKSIZE;j++)
			{
				int f_root = P[r + j] & F[p + j];
				if (!f_root)
					f_root = P[r + j];
					
				// these give the multi-state characters their final values
				// at least for now this is how those ? are treated, they get
				// final values here
				F[r + j] = f_root;
			}
#endif

		}
				
		return 0;
	}



	void reopt_source_root(cgroup_data *cd, int root, int rootHTU, int i)
	{
		st_t *prow = &cd->pstate[cd->memwidth * i];
		for (int j=0;j<BLOCKSIZE;j++)
		{
			int rv = prow[rootHTU*BLOCKSIZE+j] & 0xff; //prow[root*BLOCKSIZE+j];
			if (!rv)
			{
				rv = prow[root*BLOCKSIZE+j];
			}
			prow[root*BLOCKSIZE+j] = rv;
		}
	}
	
	// Fitch single character first pass
	int single_unordered_character_first_pass_calc_length(int *fpo, int maxnodes, int root, int rootHTU, cgroup_data *cd, int i)
	{
		int sum = 0;
		const int *bp = fpo;
		
		// offset to the right row into the submatrix table
		st_t *prow = &cd->pstate[cd->memwidth * i];
		
		#if defined(USE_SIMD)
		
			__m128i grandtot = _mm_setzero_si128();
			__m128i zero = _mm_setzero_si128();
			__m128i cmp = _mm_setzero_si128();

			int n = bp[0];
			int c1 = bp[1];
			int c2 = bp[2];
							
			while (n != -1)
			{
				__m128i a = _mm_load_si128((__m128i*)&prow[c1*BLOCKSIZE]); 
				__m128i b = _mm_load_si128((__m128i*)&prow[c2*BLOCKSIZE]);
				
				c1 = bp[3+1];
				c2 = bp[3+2];
				const int resaddr = n*BLOCKSIZE;
				
				// interleaved here, cmp starts out at z for the first round
				// and one extra at exit
				grandtot = _mm_add_epi8(grandtot, cmp);
			
				__m128i v = _mm_and_si128(a, b);
				
				cmp = _mm_cmpeq_epi8(v, zero);
				
				__m128i alt1 = _mm_and_si128(cmp, _mm_or_si128(a, b));
				__m128i alt2 = _mm_andnot_si128(cmp, v);
				__m128i res = _mm_or_si128(alt1, alt2);
				_mm_store_si128((__m128i*)&prow[resaddr], res);

				n = bp[3];
				bp += 3;
			}

			grandtot = _mm_add_epi8(grandtot, cmp);
			
			signed char tmp[BLOCKSIZE];
			_mm_storeu_si128((__m128i*)tmp, grandtot);
		
			for (int j=0;j<BLOCKSIZE;j++)
			{
				int rv = prow[rootHTU*BLOCKSIZE+j] & prow[root*BLOCKSIZE+j];
				if (!rv)
				{
					rv = prow[root*BLOCKSIZE+j];
					tmp[j]--;
				}
				prow[root*BLOCKSIZE+j] = rv;
				sum -= tmp[j] * cd->weights[i * BLOCKSIZE + j];
			}			
			
		#else

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
				
					int v = a & b;
					if (!v)
					{
						v = a | b;
						++subsum[j];
					}
					
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
		
		#endif
		return sum;
	}

	void ultranode(network::data *d, int node)
	{
		cgroup_data *cd = &d->opt->unordered;
		for (int i=0;i<cd->packunits;i++)
		{
			memset(&cd->fstate[cd->memwidth * i + BLOCKSIZE * node], 0xff, BLOCKSIZE);
		}
	}
	
	int clip_merge_dist_unordered(cgroup_data *cd, int maxnodes, int target_root, int t0, int t1, int max)
	{
		int sum = 0;
		
		// for all characters in this group
		if (t1 != network::NOT_IN_NETWORK)
		{
			target_root *= BLOCKSIZE;
			t0 *= BLOCKSIZE;
			t1 *= BLOCKSIZE;

			#if defined(USE_SIMD)
			__m128i zero = _mm_setzero_si128();
			__m128i grandtot = _mm_set1_epi32(0);
			const __m128i vk0 = _mm_set1_epi8(0);
			const __m128i vk1 = _mm_set1_epi16(1);   
			__m128i supertot;
			
			// in the hope for better caching
			if (t1 < t0)
				std::swap(t0, t1);
			
			for (int i=0;i<cd->packunits;i++)
			{
				st_t *F = &cd->fstate[cd->memwidth * i];
				__m128i ft0 = _mm_load_si128((__m128i*)&F[t0]); 
				__m128i ft1 = _mm_load_si128((__m128i*)&F[t1]);
				__m128i ftr = _mm_load_si128((__m128i*)&F[target_root]);
				__m128i wgh = _mm_load_si128((__m128i*)(&cd->weights[i * BLOCKSIZE]));
								
				__m128i rh1 = _mm_or_si128(ft0, ft1);
				__m128i tot = _mm_and_si128(ftr, rh1);
				__m128i cmp = _mm_cmpeq_epi8(tot, zero);
				__m128i subtot = _mm_and_si128(cmp, wgh);
				
				// crazy add
				__m128i vl = _mm_unpacklo_epi8(subtot, vk0);
				__m128i vh = _mm_unpackhi_epi8(subtot, vk0);
				grandtot = _mm_add_epi32(grandtot, _mm_madd_epi16(vl, vk1));
				grandtot = _mm_add_epi32(grandtot, _mm_madd_epi16(vh, vk1));
				
				if ((i & 5) == 5)
				{
					supertot = _mm_add_epi32(grandtot, _mm_srli_si128(grandtot, 8));
					supertot = _mm_add_epi32(supertot, _mm_srli_si128(supertot, 4));
					if (_mm_cvtsi128_si32(supertot) > max)
						return 100000;
				}
			}
			
			// and super sum
			supertot = _mm_add_epi32(grandtot, _mm_srli_si128(grandtot, 8));
			supertot = _mm_add_epi32(supertot, _mm_srli_si128(supertot, 4));
			sum += _mm_cvtsi128_si32(supertot);
			
			#else
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
			#endif
		}
		else
		{
			target_root *= BLOCKSIZE;
			t0 *= BLOCKSIZE;
	
			#if defined(USE_SIMD)
			__m128i zero = _mm_setzero_si128();
			__m128i grandtot = _mm_set1_epi32(0);
			const __m128i vk0 = _mm_set1_epi8(0);       
			const __m128i vk1 = _mm_set1_epi16(1);
			
			__m128i supertot;
			
			
			for (int i=0;i<cd->packunits;i++)
			{
				st_t *O = &cd->ostate[cd->memwidth * i];
				st_t *F = &cd->fstate[cd->memwidth * i];
				__m128i wgh = _mm_load_si128((__m128i*)(&cd->weights[i * BLOCKSIZE]));
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

				if ((i&3) == 0)
				{
					supertot = _mm_add_epi32(grandtot, _mm_srli_si128(grandtot, 8));
					supertot = _mm_add_epi32(supertot, _mm_srli_si128(supertot, 4));
					if (_mm_cvtsi128_si32(supertot) > max)
						return 100000;
				}
			}
			
			// here too
			grandtot = _mm_add_epi32(grandtot, _mm_srli_si128(grandtot, 8));
			grandtot = _mm_add_epi32(grandtot, _mm_srli_si128(grandtot, 4));
			sum += _mm_cvtsi128_si32(supertot);
			
			#else			
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
			#endif
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
			st_t *F = &cd->fstate[cd->memwidth * i];
			
			#if defined(USE_SIMD)
				__m128i a = _mm_load_si128((__m128i*)&F[s0]); 
				__m128i b = _mm_load_si128((__m128i*)&F[s1]);
				__m128i c = _mm_or_si128(a, b);
				_mm_store_si128((__m128i*)&F[new_node], c);
			#else
				for (int j=0;j<BLOCKSIZE;j++)
					F[new_node+j] = F[s0+j] | F[s1+j];
			#endif
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

	int optimize_for_tree(optstate *st, network::data *d, int root, bool write_final = false, bool reoptimize = false)
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


		if (!reoptimize)
		{
			// Restore taxon nodes' pstate to o-state to properly deal with multistates.
			//
			// This is not needed when reoptimizing the source tree in TBR since those values
			// will be kept (wouldn't change anyway if we actually did a first re-pass).
			for (int k=0;k<st->unordered.packunits;k++)
			{
				for (int j=0;;j++)
				{
					int p = st->first_pass_order[j];
					if (p == -1)
						break;
					
					// The P values are fiddled with during computation, so restore all for all characters that are
					// being used in this tree (we could be dealing with a subtree so don't touch other nodes data) 
					const int ofs = k * st->unordered.memwidth + p * BLOCKSIZE;
				
					if (p < d->mtx_taxons)
					{	
					#if defined(USE_SIMD)
						__m128i tmp = _mm_load_si128((__m128i*)&st->unordered.ostate[ofs]);
						_mm_store_si128((__m128i*)&st->unordered.pstate[ofs], tmp);
					#else		
						memcpy(&st->unordered.pstate[ofs], &st->unordered.ostate[ofs], BLOCKSIZE);
					#endif
					}
				}		
				
				// root isn't listed in this order list. root_htu covers case with a 2-node network where the second
				// node doesn't appear in the first_pass_order (because root has only 1 child)
				const int ofs0 = k * st->unordered.memwidth + root * BLOCKSIZE;
				const int ofs1 = k * st->unordered.memwidth + root_htu * BLOCKSIZE;
				#if defined(USE_SIMD)
					__m128i tmp0 = _mm_load_si128((__m128i*)&st->unordered.ostate[ofs0]);
					__m128i tmp1 = _mm_load_si128((__m128i*)&st->unordered.ostate[ofs1]);
					_mm_store_si128((__m128i*)&st->unordered.pstate[ofs0], tmp0);
					_mm_store_si128((__m128i*)&st->unordered.pstate[ofs1], tmp1);
				#else
					memcpy(&st->unordered.pstate[ofs0], &st->unordered.ostate[ofs0], BLOCKSIZE);
					memcpy(&st->unordered.pstate[ofs1], &st->unordered.ostate[ofs1], BLOCKSIZE);
				#endif
			}
		} // end need-to-investigate block

		DPRINT("Root=" << root << " rootHTU=" << root_htu);
		
		const int sum0 = 0;// single_unordered_character_first_pass_calc_length(st->first_pass_order, st->maxnodes, root, root_htu, &st->ordered);

		int sum1 = 0;

		if (!reoptimize)
		{
			for (int i=0;i<st->unordered.packunits;i++)
				sum1 += single_unordered_character_first_pass_calc_length(st->first_pass_order, st->maxnodes, root, root_htu, &st->unordered, i);
		}
		
		const int sum = sum0 + sum1;
		
		DPRINT("Optimized sum = " << sum0 << "+" << sum1 << "=" << sum);
		
		if (write_final)
		{
			for (int i=0;i<st->unordered.packunits;i++)
			{
				single_unordered_character_final_pass(root, st->maxnodes, d->mtx_taxons, st->net, &st->unordered, i);
			}
		}
		

		// -- lala lala --
		return sum;
	}

	character::distance_t reoptimize_final(network::data *data, network::data *ref, int root)	
	{
		network::treeify(data, root, ref->opt->net, ref->opt->first_pass_order);
		for (int i=0;i<ref->opt->unordered.packunits;i++)
			final_state_reoptimization(root, data->mtx_taxons, ref->opt->net, &data->opt->unordered, &ref->opt->unordered, i);
				
		return 0;
	}
	
	void set_weight(optstate *st, int pos, int weight)
	{
		if (pos < st->unordered.count)
		{
			st->unordered.weights[pos] = weight;
		}
		/*
		else
		{
			std::cerr << "Weight out of range, pos=" << pos << " < " << st->unordered.count << std::endl;
		}
		*/
	}

	void free(optstate *s)
	{
		delete [] s->net;
		delete s;
	}
		
		
	character::distance_t optimize(network::data *data, int root, bool write_final)
	{
		DPRINT("[optimize] - First pass run to calculate length");
		return optimize_for_tree(data->opt, data, root, write_final, 0);
	}

	void copy(optstate *target, optstate *source)
	{
		copy_cgroup(&target->ordered, &source->ordered);
		copy_cgroup(&target->unordered, &source->unordered);
	}

}
