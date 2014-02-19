#include "network.h"
#include "character.h"
#include "newick.h"
#include "optimize.h"

#include <emmintrin.h>

#include <iostream>
#include <cstring>
#include <cstdlib>

//#define DPRINT(x) { std::cout << x << std::endl; }
#define DPRINT(x) {};

namespace optimize
{
	
	namespace
	{
	
	struct fsro
	{
		st_t *F;
		st_t *P;
		st_t *FREF;
		bool *allowSkipMask;
		int taxons;
		network::node *net;
		int tw;
		int taxp[MAX_NODES];
		int taxp_count;
		int skipped, total;
		bool seq_done[MAX_CHARACTERS/BLOCKSIZE];
	};
	}
	
	void fsro_exec(fsro *s, int me, int anc)
	{
		bool *seq_done = s->seq_done;
		const int tw = s->tw;
		const int kid0 = s->net[me].c1;
		const int kid1 = s->net[me].c2;
		const int blkAncestor = tw * anc;
		const int blkMe = tw * me;
		const int blkKid0 = tw * kid0;
		const int blkKid1 = tw * kid1;
		st_t *F = s->F;
		st_t *P = s->P;
		st_t *FREF = s->FREF;
		
		int done_changes[MAX_CHARACTERS/BLOCKSIZE+1];
		int *done_changes_ptr = done_changes;
		
		for (int c=0;c<tw;c+=BLOCKSIZE)
		{
			s->total++;
			if (seq_done[c >> BLOCKSHIFT])
			{
				s->skipped++;
				continue;
			}
		
			__m128i ppMe = _mm_load_si128((__m128i*)&P[blkMe + c]);
			__m128i fa = _mm_load_si128((__m128i*)&F[blkAncestor + c]);
			__m128i pp0 = _mm_load_si128((__m128i*)&P[blkKid0 + c]);
			__m128i pp1 = _mm_load_si128((__m128i*)&P[blkKid1 + c]);
			__m128i fme = _mm_and_si128(ppMe, fa);
			__m128i fref = _mm_load_si128((__m128i*)&FREF[blkMe + c]);
			
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
		
			_mm_store_si128((__m128i*)&F[blkMe + c], totalResult);
			
			__m128i vcmp = (__m128i)_mm_cmpeq_ps((__m128)fref, (__m128)totalResult);
			
			if (_mm_movemask_epi8(vcmp) == 0xffff)
			{
				if (!s->allowSkipMask || s->allowSkipMask[me])
				{
					unsigned int idx = (c >> BLOCKSHIFT);
					if (!seq_done[idx])
					{
						seq_done[idx] = true;
						*done_changes_ptr++ = idx;
					}
				}
			}
		}
		
		if (kid0 >= s->taxons)
		{
			fsro_exec(s, kid0, me);
		}
		else if (kid0 >= 0)
		{
			s->taxp[s->taxp_count] = kid0;
			s->taxp[s->taxp_count+1] = me;
			s->taxp_count += 2;
		}
		
		if (kid1 >= s->taxons)
		{
			fsro_exec(s, kid1, me);
		}
		else if (kid1 >= 0)
		{
			s->taxp[s->taxp_count] = kid1;
			s->taxp[s->taxp_count+1] = me;
			s->taxp_count += 2;
		}
		
		// undo the done mask
		while (done_changes_ptr > done_changes)
		{
			done_changes_ptr--;
			seq_done[*done_changes_ptr] = false;
		}
	}

	// Fitch single character final pass
	void final_state_reoptimization(int root, int taxons, network::node *net, cgroup_data *cd, cgroup_data *ref, bool is_tree_root=true, bool *allowSkipMask=0)
	{
		
		fsro fs;
		
		fs.tw = cd->taxonwidth;
		fs.F = cd->fstate;
		fs.P = cd->pstate;
		fs.FREF = ref->fstate;
		fs.net = net;
		fs.skipped = 0;
		fs.total = 0;
		fs.taxons = taxons;
		fs.allowSkipMask = allowSkipMask;
		memset(fs.seq_done, 0x00, sizeof(fs.seq_done));

		const int tw = cd->taxonwidth;
		
		if (is_tree_root)
		{
			memcpy(&fs.F[tw * root], &fs.P[tw * root], tw);
			
			fs.taxp[0] = root;
			fs.taxp_count = 2;
			
			if (net[root].c0 >= 0)
			{
				fsro_exec(&fs, net[root].c0, root);
				fs.taxp[1] = net[root].c0;
			}
			else if (net[root].c1 >= 0)
			{
				fsro_exec(&fs, net[root].c1, root);
				fs.taxp[1] = net[root].c1;
			}
			else if (net[root].c2 >= 0)
			{
				fsro_exec(&fs, net[root].c2, root);
				fs.taxp[1] = net[root].c2;
			}
		}
		else
		{
			memcpy(&fs.F[tw * root], &fs.FREF[tw * root], tw);
	
			// need to put things in right queue
			int kids[2] = {net[root].c1, net[root].c2};
			
			fs.taxp_count = 0;
			
			for (int k=0;k<2;k++)
			{
				if (kids[k] >= taxons)
				{
					fsro_exec(&fs, kids[k], root);
				}
				else if (kids[k] >= 0)
				{
					fs.taxp[fs.taxp_count] = kids[k];
					fs.taxp[fs.taxp_count+1] = root;
					fs.taxp_count += 2;
				}
			}
		}
		
		// Need to treat terminal nodes specially when they can have ? in them
		for (int A=0;A<fs.taxp_count;A+=2)
		{
			const int r = tw * fs.taxp[A];
			const int p = tw * fs.taxp[A+1];
			//std::cout << "w[" << fs.taxp[A] << "] from " << fs.taxp[A+1] << std::endl;
			for (int c=0;c<tw;c+=BLOCKSIZE)
			{
				__m128i zero = _mm_setzero_si128();
				__m128i PR = _mm_load_si128((__m128i*)&fs.P[r + c]);
				__m128i FP = _mm_load_si128((__m128i*)&fs.F[p + c]);
				__m128i f_root = _mm_and_si128(PR, FP);
				__m128i cmp = _mm_cmpeq_epi8(f_root, zero);
				__m128i totalResult = _mm_or_si128(_mm_andnot_si128(cmp, f_root),
									     _mm_and_si128(cmp, PR));
				_mm_store_si128((__m128i*)&fs.F[r + c], totalResult);
			}
		}
	}

	// Fitch single character first pass
	int single_unordered_character_first_pass_reoptimize(int *fpo, int maxnodes, int root, int rootHTU, cgroup_data *cd, cgroup_data *ref)
	{
		const int *bp = fpo;
		
		// offset to the right row into the submatrix table
		st_t *prow = cd->pstate;
		st_t *refrow = ref->pstate;
		
		int tw = cd->taxonwidth;
		int lastChanged = -1;

		int n = bp[0];
		int c1 = bp[1];
		int c2 = bp[2];

		__m128i zero = _mm_setzero_si128();
							
		while (n != -1)
		{
			for (int cp=0;cp<tw;cp+=BLOCKSIZE)
			{
				const int resaddr = tw * n + cp;
				__m128i a = _mm_load_si128((__m128i*)&prow[c1*tw+cp]);
				__m128i b = _mm_load_si128((__m128i*)&prow[c2*tw+cp]);
				__m128i ref = _mm_load_si128((__m128i*)&refrow[resaddr]);
				__m128i v = _mm_and_si128(a, b);
				__m128i cmp = _mm_cmpeq_epi8(v, zero);
				__m128i alt1 = _mm_and_si128(cmp, _mm_or_si128(a, b));
				__m128i alt2 = _mm_andnot_si128(cmp, v);
				__m128i res = _mm_or_si128(alt1, alt2);
				_mm_store_si128((__m128i*)&prow[resaddr], res);
				
				// results now agree
				__m128i vcmp = (__m128i)_mm_cmpeq_ps((__m128)ref, (__m128)res);
				if (_mm_movemask_epi8(vcmp) != 0xffff || lastChanged == -1)
					lastChanged = n;
			}
			
			n = bp[3];
			c1 = bp[3+1];
			c2 = bp[3+2];
			bp += 3;
		}

		for (int j=0;j<tw;j++)
		{
			int rv = prow[rootHTU*tw+j] & prow[root*tw+j];
			if (!rv)
				rv = prow[root*tw+j];
			prow[root*tw+j] = rv;
		}
	
		return lastChanged;
	}

	void tbr_target_reoptimization(network::data *data, network::data *ref, int calcroot, int start)
	{
	
		int fpo[1024];
		int *fpo_out = fpo;
		
		optstate *st = data->opt;		
		network::node *net = st->net;
		network::treeify(data, calcroot, net, fpo);

		
		if (start < data->mtx_taxons)
		{
			if (net[start].c0 != -1)
				start = net[start].c0;
			else if (net[start].c1 != -1)
				start = net[start].c1;
			else
				start = net[start].c2;
		}

		bool skipMask[1024];
		memset(skipMask, 0x1, data->allocnodes);

		// construct traverse list from clip point all the way to the tree for 1st pass reopt
		int cur = start;
		while (cur >= 0 && cur != calcroot)
		{
			skipMask[cur] = false;
			fpo_out[0] = cur;
			fpo_out[1] = net[cur].c1;
			fpo_out[2] = net[cur].c2;
			fpo_out += 3;
			cur = net[cur].c0;		
		}
		
		skipMask[calcroot] = false;

		// ---
		*fpo_out = -1;

		//
		const int root = calcroot;
		int root_htu = net[root].c1;
		if (root_htu < 0)
			root_htu = net[root].c2;

		reset_taxons_states(&st->unordered, fpo, data->mtx_taxons, root, root_htu);

		// go up to the root
		int topmost_change = single_unordered_character_first_pass_reoptimize(fpo, st->maxnodes, root, root_htu, &st->unordered, &ref->opt->unordered);

		if (topmost_change == -1)
			topmost_change = start;
		
		if (net[topmost_change].c0 != -1)
			topmost_change = net[topmost_change].c0;

		skipMask[topmost_change] = false;

		if (net[topmost_change].c0 != -1)
			topmost_change = net[topmost_change].c0;

		skipMask[topmost_change] = false;
		
		final_state_reoptimization(topmost_change, data->mtx_taxons, st->net, &data->opt->unordered, &ref->opt->unordered, topmost_change == calcroot, skipMask);
	}
	
	void tbr_source_reoptimization(network::data *data, network::data *ref, int root)
	{
		network::treeify(data, root, ref->opt->net, ref->opt->first_pass_order);
		final_state_reoptimization(root, data->mtx_taxons, ref->opt->net, &data->opt->unordered, &ref->opt->unordered);
	}
}
