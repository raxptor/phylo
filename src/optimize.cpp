#include "network.h"
#include "character.h"
#include "newick.h"

#include <iostream>
#include <cstring>

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
				
				out->pstate[dataofs] = out->fstate[dataofs] = mask(source->submatrix[mtxofs]);
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
	
	struct unordered_scoring 
	{
		static inline void mod_score(int unionmask, int *out)
		{
			++(*out);
		}
	};
	
	struct ordered_scoring
	{
		static inline void mod_score(int unionmask, int *out)
		{
			(*out) += mindist[unionmask];
		}
	};
	
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
				std::cout << to_char(cd->fstate[i * maxnodes + t]);
			}
			std::cout << " p => ";
			for (unsigned int i=0;i<cd->count;i++)
			{
				std::cout << to_char(cd->pstate[i * maxnodes + t]);
			}
			std::cout << std::endl;
		}
	}
	
	template<typename SCORING_FN>
	int slow_single_first_pass(int *fpo, int maxnodes, int root, int rootHTU, cgroup_data *cd)
	{
		int sum = 0;
		
		// for all characters in this group
		for (int i=0;i<cd->count;i++)
		{
			const int *bp = fpo;

			// offset to the right row into the submatrix table
			st_t *prow = &cd->pstate[maxnodes * i];

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
					SCORING_FN::mod_score(v, &sum);
				}
				
				DPRINT("writing pstate[" << n << "] [char:" << i << "] = " << v << " sum=" << sum);
				prow[n] = v;
			}
			
			int rv = prow[rootHTU] & prow[root];
			if (!rv)
			{
				rv = prow[rootHTU] | prow[root];
				SCORING_FN::mod_score(prow[rootHTU] | prow[root], &sum);
				DPRINT("Diff at root, scor=" << sum);
			}
			
			DPRINT("offsetroot = " << (cd->count * i + root));
			DPRINT("root@" << root << " => " << rv);
			prow[root] = rv;
		}
		
		return sum;	
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
	
	int optimize_for_tree(optstate *st, network::data *d, int root)
	{
		
		st->root = root;
		network::treeify(d, root, st->net, st->first_pass_order);
		
		int rootHTU = st->net[root].c1;
		if (rootHTU < 0)
			rootHTU = st->net[root].c2;
			
		DPRINT("Root=" << root << " rootHTU=" << rootHTU);

		int sum0, sum1;
		sum0 = slow_single_first_pass<unordered_scoring>(st->first_pass_order, st->maxnodes, root, rootHTU, &st->ordered);
		sum1 = slow_single_first_pass<unordered_scoring>(st->first_pass_order, st->maxnodes, root, rootHTU, &st->unordered);
		int sum = sum0 + sum1;

//		print_state(st, st->maxnodes, d->mtx_taxons);

		DPRINT("Optimized sum = " << sum0 << "+" << sum1 << "=" << sum);
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


	/*
	void push_down(optstate *s, int w_, char value_)
	{
		struct que {
			int w;
			char value;
		};
		
		que next[1024];
		int queue = 0;
		
		next[0].w = w_;
		next[0].value = value_;

		while (queue >= 0)
		{
			const int w = next[queue].w;
			const int value = next[queue].value;
			if (w >= s->taxons)
			{
				const char sv = is_single[(unsigned char)s->bmp[w]];
				character::state_t write;
				if (sv != -1)
				{
					// definite value
					DPRINT("writing single value " << (int)sv << " because mask " << (int)s->bmp[w]);
					write = sv;
					s->sum += dist(value, sv);
				}
				else
				{
					char mv = mask(value);
					if ((mv & s->bmp[w]))
					{
						DPRINT("writing single value " << (int)value << " because it matched the mask " << (int)s->bmp[w]);
						write = value;
					}
					else
					{
						for (int i=0;i<32;i++)
						{
							if (mask(i) & s->bmp[w])
							{
								write = i;
								break;
							}
						}
						DPRINT("writing " << write << " because writemask");
					}
				}
				
				if (value != character::UNKNOWN_CHAR_VALUE && value != write)
				{
					DPRINT("step from " << int(value) << " to " << int(write));
					s->sum += dist(value, write);
				}
				
				// remove one, add one
				next[queue].w = s->net[w].c1;
				next[queue].value = write;
				next[queue+1].w = s->net[w].c2;
				next[queue+1].value = write;
				s->cval[w] = write;
				queue++;
				continue;
			}
				
			if (w == network::NOT_IN_NETWORK)
			{
				queue--;
				continue;
			}
		
			if (s->cval[w] == character::UNKNOWN_CHAR_VALUE)
			{
				queue--;
				continue;
			}
				
			if (value != s->cval[w])
			{
				DPRINT("step from " << int(value) << " to " << int(s->cval[w]));
				s->sum += dist(value, s->cval[w]);
			}
			
			queue--;
		}
	}
	*/
	
	character::distance_t optimize(network::data *data, bool all_chars)
	{
		DPRINT("[optimize] - Full optimization run.");
		return optimize_for_tree(data->opt, data, 0);
		
		
	}
}
