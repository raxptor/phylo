#include "network.h"
#include "character.h"
#include <iostream>


//#define DPRINT(x) { std::cout << x << std::endl; }
#define DPRINT(x) {};

namespace optimize
{

	struct state
	{
		int root;
		int character;
		char *bmp;
		character::state_t **cdata;
		network::node *net;
		int taxons;
	};
	
	char is_single[256];

	char mask(int val)
	{
		return 1 << val;
	}

	void down(network::data *data, int where, int depth)
	{
		if (where == network::NOT_IN_NETWORK)
			return;
			
		std::cout << where << " " << depth << std::endl;
		down(data, data->network[where].c1, depth + 1);
		down(data, data->network[where].c2, depth + 1);
	}
	
	char collect_two(state *s, int w)
	{
		if (w == network::NOT_IN_NETWORK)
			return 0;
		
		if (w < s->taxons && w != s->root)
		{
			return mask(s->cdata[w][s->character]);
		}
			
		char a = collect_two(s, s->net[w].c1);
		char b = collect_two(s, s->net[w].c2);
		char tmp;
		
		// if intersection exists, use that, otherwise union
		tmp = a & b;

		if (!tmp)
			tmp = a | b;

		s->bmp[w] = tmp;
		return tmp;		
	}
	
	void push_down(state *s, int w, char value)
	{
		if (w == network::NOT_IN_NETWORK)
			return;
		if (w < s->taxons)
			return;
		
		char sv = is_single[s->bmp[w]];
		if (sv != -1)
		{
			// definite value
			s->cdata[w][s->character] = sv;
		}
		else
		{
			char mv = mask(value);
			if (mv & s->bmp[w])
			{
				s->cdata[w][s->character] = value;
			}
			else
			{
				for (int i=0;i<32;i++)
					if (mask(i) & s->bmp[w])
					{
						s->cdata[w][s->character] = i;
						if (i == 31)
							std::cout << "Found nothing!" << std::endl;
						break;
					}
			}
		}
		push_down(s, s->net[w].c1, s->cdata[w][s->character]);
		push_down(s, s->net[w].c2, s->cdata[w][s->character]);
	}

	void optimize(network::data *data)
	{
		DPRINT("Treeifying network...");
		network::treeify(data, 0);
		
		state s;
		s.bmp = new char[data->allocnodes];
		
		for (int i=0;i<256;i++)
			is_single[i] = -1;
		for (int i=0;i<8;i++)
			is_single[0 << i] = i;
			
		s.root = 0;
		s.net = data->network;
		s.taxons = data->mtx_taxons;
		s.cdata = data->characters;

		for (int i=0;i<data->mtx_characters;i++)
		{
			s.character = i;
			memset(s.bmp, 0x00, data->allocnodes);
			collect_two(&s, s.root);
			s.bmp[s.root] = mask(s.cdata[s.root][i]);
			push_down(&s, s.net[s.root].c1, s.cdata[s.root][i]);
			push_down(&s, s.net[s.root].c2, s.cdata[s.root][i]);
		}
		
		delete[] s.bmp;
//		std::cout << "Previous distance:" << data->dist << std::endl;
		network::recompute_dist(data);
//		std::cout << "New distance:" << data->dist << std::endl;
	}
}
