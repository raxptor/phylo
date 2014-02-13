#include "network.h"
#include "character.h"
#include "newick.h"

#include <iostream>
#include <cstring>

//#define DPRINT(x) { std::cout << x << std::endl; }
#define DPRINT(x) {};

namespace optimize
{
	// character
	struct ocstate
	{
		character::distance_t sum;
	};
	
	struct optstate
	{
		ocstate *chars;
		int characters;
		
		// work memory
		int root;
		int character;
		char *bmp;
		network::node *net;
		int sum;
		int totsum;
		int taxons;
		int netsize;
		int bottomup[1024];
		
		int pstate[1024];		
		
		int cval[1024];
	};

	void copy(optstate *target, optstate *source)
	{
		memcpy(target->chars, source->chars, sizeof(ocstate) * source->characters);
	}
	
	optstate* create(int numchars, int allocnodes)
	{
		optstate *st = new optstate();
		st->chars = new ocstate[numchars];
		st->bmp = new char[allocnodes];
		st->net = new network::node[allocnodes];
		st->characters = numchars;
		
		for (int i=0;i<numchars;i++)
			st->chars[i].sum = 0;
		
		return st;
	}
	
	void free(optstate *s)
	{
		delete [] s->net;
		delete [] s->bmp;
		delete [] s->chars;
		delete s;
	}
	
	int is_single[256];
	int mindist[256];

	void init()
	{
		for (int i=0;i<256;i++)
			is_single[i] = -1;
			
		for (int i=0;i<32;i++)
			is_single[0 << i] = i;
			
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

	inline char mask(int val)
	{
		if (val == character::UNKNOWN_CHAR_VALUE)
			return 0xff;
			
		return 1 << val;
	}
	
	// these are all valid values
	inline character::distance_t dist(character::state_t a, character::state_t b)
	{
		if (a > b) 
			return a - b; 
		return b - a;
	}
	
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
	
	character::distance_t optimize(network::data *data, bool all_chars)
	{
		DPRINT("Treeifying network...");
		optstate & s = *(data->opt);
		
		// unnecessary but why not
		
		s.root = 0;
		s.taxons = data->mtx_taxons;
		s.totsum = 0;

		network::treeify(data, 0, s.net, s.bottomup);
//		newick::print(data);
		
		s.netsize = 1;
		for (int i=0;i<2048;i++)
		{
			if (s.bottomup[i] < 0) 
			{
				s.netsize = i;
				break;
			}
		}
		
		if (s.netsize % 3)
		{
			std::cerr << "Visit order is not %3=0, it is " << s.netsize << "!" << std::endl;
			exit(-2);
		}
		
		DPRINT("Netsize = " << s.netsize);		
		
		int sum = 0;
		
		for (int c=0;c<data->mtx_characters;c++)
		{
			for (int i=0;i<data->mtx_taxons;i++)
				s.pstate[i] = mask(data->characters[i][c]);
			
			for (int i=0;i<s.netsize;i+=3)
			{
				const int n = s.bottomup[i];
				const int c1 = s.bottomup[i+1];
				const int c2 = s.bottomup[i+2];

				const int a = s.pstate[c1];
				const int b = s.pstate[c2];
				
				int v = a & b;
				if (!v)
				{
					v = a | b;
					DPRINT("Visiting " << n << " c1=" << c1 << " c2=" << c2 << " union " << (int)a << " | " << (int)b << " score=" << mindist[v]);
					sum += mindist[v];
				}
				
				DPRINT("writing pstate[" << n << "] [" << i << "/" << s.netsize << "] = " << v);
				s.pstate[n] = v;
				s.bmp[n] = v;
			}
			
			int rootKid = s.net[s.root].c1;
			if (rootKid < 0)
				rootKid = s.net[s.root].c2;
				
			int rv = s.pstate[rootKid] & s.pstate[s.root];
			if (!rv)
			{
				DPRINT(" root adds some too " << (int)(s.pstate[rootKid] | s.pstate[s.root]) << " " << (int)rootKid << " " << mindist[s.pstate[rootKid] | s.pstate[s.root]]);
				DPRINT(" mindislookup = " << (int)(s.pstate[rootKid]|s.pstate[s.root]));
				DPRINT(" pstate root = " << s.pstate[s.root]);
				DPRINT(" pstate kid = " << s.pstate[rootKid] << " slot=" << rootKid);
				sum += mindist[s.pstate[rootKid] | s.pstate[s.root]];
			}
		}
		
		DPRINT("Optimized sum = " << sum);
		// -- lala lala --
		return sum;
	}
}
