#include "network.h"
#include "character.h"
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
		int taxons;
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

	void init()
	{
		for (int i=0;i<256;i++)
			is_single[i] = -1;
		for (int i=0;i<32;i++)
			is_single[0 << i] = i;
		is_single[0 << character::UNKNOWN_CHAR_VALUE] = -1;
	}

	inline char mask(int val)
	{
		if (val == character::UNKNOWN_CHAR_VALUE)
			return 0x0;
			
		return 1 << val;
	}
	
	// these are all valid values
	character::distance_t dist(character::state_t a, character::state_t b)
	{
		if (a > b) 
			return a - b; 
		return b - a;
	}

	void down(network::data *data, int where, int depth)
	{
		if (where == network::NOT_IN_NETWORK)
			return;
			
		std::cout << where << " " << depth << std::endl;
		down(data, data->network[where].c1, depth + 1);
		down(data, data->network[where].c2, depth + 1);
	}
	
	char collect_two(optstate *s, int w)
	{
		const network::node *net = s->net;
		const int c1 = net[w].c1;
		const int c2 = net[w].c2;
		
		if (w < s->taxons && w != s->root)
		{
			DPRINT("mask[" << w << "] is " << (int)mask(s->cval[w]));
			return mask(s->cval[w]);
		}
		
		char a=0, b=0;
		if (c1 != network::NOT_IN_NETWORK)
			a = collect_two(s, c1);
		if (c2 != network::NOT_IN_NETWORK)
			b = collect_two(s, c2);
			
		// if intersection exists, use that, otherwise union
		char tmp = a & b;

		if (!tmp)
			tmp = a | b;

		DPRINT("mask[" << w << "] is combination " << tmp);

		s->bmp[w] = tmp;
		return tmp;		
	}
	
	void push_down(optstate *s, int w_, char value_)
	{
		DPRINT("traverse [" << w << "] val=" << (int)value);
		
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
				const char sv = is_single[s->bmp[w]];
				character::state_t write;
				if (sv != -1)
				{
					// definite value
					DPRINT("writing single value " << (int)sv << " because mask " << (int)s->bmp[w]);
					write = sv;
					if (sv != value)
					{
						DPRINT("step from " << (int)sv << " to " << (int)value);
						s->sum += dist(value, sv);
					}
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
	
	void run_character(network::data *data, optstate *s, int i)
	{
		s->sum = 0;
		s->character = i;
		
		collect_two(s, s->root);
		
		s->bmp[s->root] = mask(s->cval[s->root]);
		push_down(s, s->net[s->root].c1, s->cval[s->root]);
		push_down(s, s->net[s->root].c2, s->cval[s->root]);
		data->opt->chars[i].sum = s->sum;
	}
	
	character::distance_t optimize(network::data *data, bool all_chars)
	{
		DPRINT("Treeifying network...");
		optstate & s = *(data->opt);
		
		// unnecessary but why not
		s.root = 0;
		s.taxons = data->mtx_taxons;

		network::treeify(data, 0, s.net);
		
		int count = 0;
		for (int i=0;i<data->mtx_characters;i++)
		{
			// todo fix charactetrs per row
			for (int j=0;j<data->allocnodes;j++)
				s.cval[j] = data->characters[j][i];
		
			run_character(data, &s, i);

			for (int j=0;j<data->allocnodes;j++)
				data->characters[j][i] = s.cval[j];
		}
		
		DPRINT("Recomputed " << count << " characters");
		
		character::distance_t sum = 0;
		for (int i=0;i<data->mtx_characters;i++)
			sum += data->opt->chars[i].sum;

/*			
		int sum2 = network::distance_by_edges(data);	
		if (sum != sum2)
		{
			network::print_characters(data);
			std::cerr << "Fucked up the count again " << sum2 << " but i say " << sum << std::endl;
			exit(1);
		}
*/
		DPRINT("  [opt] - prevdist=" << data->dist << " newdist=" << sum);
		return sum;
	}
}
