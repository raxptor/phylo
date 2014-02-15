#include "network.h"
#include "character.h"
#include "optimize.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>

// Debug stuffs
//#define NETWORKDEBUG
//#define NETWORK_CHECKS

#if defined(NETWORKDEBUG)
	#define DPRINT(x) { std::cout << x << std::endl; }
#else
	#define DPRINT(x) { }
#endif

namespace 
{
	// it prints negative numbers as x
	void itoa_hex(int num, char *where)
	{
		static const char *hex = "0123456789abcdef";
		if (num < 0)
		{
			*where++ = 'x';
		}
		else if (num < 0x10)
		{
			*where++ = hex[num];
		}
		else if (num < 0x100)
		{
			where[0] = hex[(num >> 4) & 0xf];
			where[1] = hex[num & 0xf];
			where += 2;
		}
		else if (num < 0x1000)
		{
			where[0] = hex[(num >> 8) & 0xf];
			where[0] = hex[(num >> 4) & 0xf];
			where[1] = hex[num & 0xf];
			where += 3;
		}
		else if (num < 0x10000)
		{
			where[0] = hex[(num >> 12) & 0xf];
			where[1] = hex[(num >> 8) & 0xf];
			where[2] = hex[(num >> 4) & 0xf];
			where[3] = hex[num & 0xf];
			where += 4;
		}
		else
		{
			std::cerr << "Fix itoa_hex. " << num << " too large" << std::endl;
		}
		*where = 0;
	}
}

namespace network
{
	data* alloc(matrix::data *mtx)
	{
		data *d = new data();
		d->matrix = mtx;
		d->mtx_taxons = mtx->taxons;
		d->mtx_characters = mtx->characters;
		d->allocnodes = 2 * mtx->taxons - 2;
		d->network = new node[d->allocnodes];

		// free list
		d->freecount = d->allocnodes - mtx->taxons;
		d->freelist = new idx_t[d->freecount];
		for (unsigned int i=0;i<d->freecount;i++)
			d->freelist[i] = i + mtx->taxons;

		d->opt = optimize::create(d);
		return d;
	}
	
	// 
	void copy(data *target, data *source)
	{
		memcpy(target->network, source->network, sizeof(node) * source->allocnodes);
		optimize::copy(target->opt, source->opt);
		memcpy(target->freelist, source->freelist, sizeof(idx_t) * target->freecount);
		target->freecount = source->freecount;
	}

	void free(data *d)
	{
		delete [] d->network;
		optimize::free(d->opt);
		delete d;
	}

	void init(network::data *d, idx_t taxon0, idx_t taxon1)
	{
		for (unsigned int i=0;i<d->allocnodes;i++)
		{
			d->network[i].c0 = NOT_IN_NETWORK;
			d->network[i].c1 = NOT_IN_NETWORK;
			d->network[i].c2 = NOT_IN_NETWORK;
		}

		d->freecount = d->allocnodes - d->matrix->taxons;
		for (unsigned int i=0;i<d->freecount;i++)
			d->freelist[i] = i + d->matrix->taxons;
	
		d->network[taxon0].c0 = taxon1;
		d->network[taxon1].c0 = taxon0;
		
		DPRINT("Network starts with " << taxon0 << " and " << taxon1 << std::endl);
	}

	void check(network::data *d)
	{
#if defined(NETWORK_CHECKS)
#endif
	}

	// leaves c2 on newwmiddle unconnected.	
	void edge_split(node *net, idx_t n0, idx_t n1, idx_t newmiddle)
	{
		node * const N0 = &net[n0];
		node * const N1 = &net[n1];

		if (N0->c0 == n1)
			N0->c0 = newmiddle;
		else if (N0->c1 == n1)
			N0->c1 = newmiddle;
		else if (N0->c2 == n1)
			N0->c2 = newmiddle;
			
		if (N1->c0 == n0)
			N1->c0 = newmiddle;
		else if (N1->c1 == n0)
			N1->c1 = newmiddle;
		else if (N1->c2 == n0)
			N1->c2 = newmiddle;
			
		net[newmiddle].c0 = n0;
		net[newmiddle].c1 = n1;
	}
	
	void edge_merge(node *net, idx_t n0, idx_t n1, idx_t middle)
	{
		node * const R0 = &net[n0];
		node * const R1 = &net[n1];
		
		// r0 & o2 to merge
		if (R0->c0 == middle)
			R0->c0 = n1;
		else if (R0->c1 == middle)
			R0->c1 = n1;
		else if (R0->c2 == middle)
			R0->c2 = n1;

		if (R1->c0 == middle)
			R1->c0 = n0;
		else if (R1->c1 == middle)
			R1->c1 = n0;
		else if (R1->c2 == middle)
			R1->c2 = n0;
	}
	
	idx_t insert(data *d, idx_t n0, idx_t n1, idx_t which)
	{
		DPRINT(" Insert " << which << " at " << n0 << "-" << n1);

		node * const net = d->network;
		
		const idx_t n = node_alloc(d);

		// insert middle		
		edge_split(d->network, n0, n1, n);

		net[n].c2 = which;
		
		if (which < d->mtx_taxons)
		{
			// this will be leaf node
			net[which].c0 = n;
			net[which].c1 = NOT_IN_NETWORK;
			net[which].c2 = NOT_IN_NETWORK;
		}
		else
		{
			DPRINT("     (inserting subtree)");
			
			// must have all the remaining connection at c1, c2
			net[which].c0 = n;
		}

		return n;
	}
	
	void the_two_others(node *network, idx_t where, idx_t which, idx_t *r0, idx_t *r1)
	{
		if (network[where].c0 == which)
		{
			*r0 = network[where].c1;
			*r1 = network[where].c2;
		}
		else if (network[where].c1 == which)
		{
			*r0 = network[where].c0;
			*r1 = network[where].c2;
		}
		else if (network[where].c2 == which)
		{
			*r0 = network[where].c0;
			*r1 = network[where].c1;
		}
		else
		{
			std::cerr << "the-two-others error: node " << which << " is no neighbour of " << where << std::endl;
		}
	}

	void disconnect(data *d, idx_t which)
	{
		node * network = d->network;

		const idx_t n = network[which].c0;
		DPRINT("Disconnecting " << which << " with 'parent' " << n);
				
		idx_t r0, r1;
		the_two_others(network, n, which, &r0, &r1);
		
		DPRINT(" -> merging from " << r0 << "-" << r1 << " with " << n << " in the middle");
		edge_merge(d->network, r0, r1, n);

		// which must a terminal node
		node_free(d, n);

#if defined(NETWORK_CHECKS)
		check(d);
#endif
	}
	
	idx_t node_alloc(data *d)
	{
		if (!d->freecount)
		{
			std::cerr << "[HELP] extra alloc failed" << std::endl;
			exit(-2);
		}
		
		return d->freelist[--d->freecount];
	}
	
	void node_free(data *d, idx_t where)
	{
		if (where < d->matrix->taxons)
		{
			std::cerr << "Trying to free original taxon" << std::endl;
		}
		
		d->freelist[d->freecount++] = where;
	}
	
	void print_characters(data *d)
	{
		std::cout << "print_characters FIXME" << std::endl;
		for (int i=0;i<d->allocnodes;i++)
		{
			// std::cout << i << " " << character::to_string(d->matrix->taxonbase[i], d->mtx_characters) << std::endl;
		}
	}

	void trace_edgelist(data *d, idx_t start, edgelist * out)
	{
		node *net = d->network;
		idx_t *outptr = out->pairs;
		int queue = 0;
		
		idx_t toexplore[MAX_NODES];
		idx_t source[MAX_NODES];

		toexplore[0] = start;
		source[0] = NOT_IN_NETWORK;
		
		while (queue >= 0)
		{
			const idx_t cur = toexplore[queue];
			const idx_t src = source[queue];
			--queue;

			const idx_t c[3] = { net[cur].c0, net[cur].c1, net[cur].c2 };
			
			for (int i=0;i<3;i++)
			{
				if (c[i] != src && c[i] != NOT_IN_NETWORK)
				{
					outptr[0] = cur;
					outptr[1] = c[i]; 
					outptr += 2;
					++queue;
					toexplore[queue] = c[i];
					source[queue] = cur;
				}
			}
		}
		
		out->count = (outptr - out->pairs);
	}

	void make_traverse_order(network::data *data, idx_t root, idx_t *bottomup, idx_t *root_htu)
	{
		// Assumes the root node is ate least connected to 1 htu, else it breaks
		idx_t toexplore[MAX_NODES];
		idx_t source[MAX_NODES];
		node *net = data->network;

		// put [root_htu], [root] on the above stack
		source[0] = root;
		source[1] = NOT_IN_NETWORK;
		toexplore[1] = root;

		// write root htu
		if (net[root].c0 != NOT_IN_NETWORK)
			toexplore[0] = *root_htu = net[root].c0;
		else if (net[root].c1 != NOT_IN_NETWORK)
			toexplore[0] = *root_htu = net[root].c1;
		else 
			toexplore[0] = *root_htu = net[root].c2;

		int queue = 1;
		int tmpOrder[1024];
		int tmpOrderOut = 0;
		while (queue >= 0)
		{
			const idx_t cur = toexplore[queue];
			const idx_t src = source[queue];
			--queue;
			
			// only HTU & root
			if (cur >= data->mtx_taxons)
			{
				tmpOrder[tmpOrderOut++] = cur;
				if (net[cur].c0 != src)
				{
					++queue;
					source[queue] = cur;
					toexplore[queue] = tmpOrder[tmpOrderOut] = net[cur].c0;
					tmpOrderOut++;
				}
				if (net[cur].c1 != src)
				{
					++queue;
					source[queue] = cur;
					toexplore[queue] = tmpOrder[tmpOrderOut] = net[cur].c1;
					tmpOrderOut++;
				}
				if (net[cur].c2 != src)
				{
					++queue;
					source[queue] = cur;
					toexplore[queue] = tmpOrder[tmpOrderOut] = net[cur].c2;
					tmpOrderOut++;
				}
			}
		}

		// reverse it		
		for (int i=0;i<tmpOrderOut;i+=3)
		{
			int up = tmpOrderOut - 3 - i;
			bottomup[up] = tmpOrder[i];
			bottomup[up+1] = tmpOrder[i+1];
			bottomup[up+2] = tmpOrder[i+2];
		}
		bottomup[tmpOrderOut] = -1;
	}

	void to_string(network::data *data, char *buffer, unsigned int bufsize)
	{
		/*
		if (data->freecount != 0)
		{
			strcpy(buffer, "incomplete net");
			return;
		}*/
		
		char *ptr = buffer;
		*ptr = 0;
		for (int i=0;i<data->allocnodes;i++)
		{
			if (ptr - buffer > (bufsize - 32))
			{
				*ptr = 0;
				std::cerr << "ERROR: Too small character buffer to to_string" << std::endl;
				exit(-2);
			}
			char buf[3][32];
			itoa_hex(data->network[i].c0, buf[0]);
			itoa_hex(data->network[i].c1, buf[1]);
			itoa_hex(data->network[i].c2, buf[2]);
			strcpy(ptr, buf[0]); 
			while (*ptr) ++ptr;
			*ptr++ = '-';
			strcpy(ptr, buf[1]);
			while (*ptr) ++ptr;
			*ptr++ = '-';
			strcpy(ptr, buf[2]);
			while (*ptr) ++ptr;
			*ptr++ = '|';
		}
		*ptr = 0;
	}
	
	void treeify(network::data *data, idx_t root, node *out, idx_t *bottomup)
	{
		node *net = data->network;
		
		if (!out)
			out = net;
			
		int queue = 0;
		
		idx_t toexplore[MAX_NODES];
		idx_t source[MAX_NODES];

		toexplore[0] = root;
		source[0] = NOT_IN_NETWORK;
		
		int tmpOrder[1024];
		int tmpOrderOut = 0;
		
		while (queue >= 0)
		{
			const idx_t cur = toexplore[queue];
			const idx_t src = source[queue];
			--queue;
			
			// only HTU & root
			if (cur >= data->mtx_taxons)
			{
				tmpOrder[tmpOrderOut++] = cur;
				if (net[cur].c0 != src)
					tmpOrder[tmpOrderOut++] = net[cur].c0;		
				if (net[cur].c1 != src)
					tmpOrder[tmpOrderOut++] = net[cur].c1;
				if (net[cur].c2 != src)
					tmpOrder[tmpOrderOut++] = net[cur].c2;
			}
			
			// always make c0 point upwards
			out[cur] = net[cur];
			
			if (out[cur].c1 == src)
				std::swap(out[cur].c0, out[cur].c1);
			else if (out[cur].c2 == src)
				std::swap(out[cur].c0, out[cur].c2);
				
			const idx_t c[3] = { net[cur].c0, net[cur].c1, net[cur].c2 };
			for (int i=0;i<3;i++)
			{
				if (c[i] != src && c[i] != NOT_IN_NETWORK)
				{
					++queue;
					toexplore[queue] = c[i];
					source[queue] = cur;
				}
			}
		}
		
		for (int i=0;i<tmpOrderOut;i+=3)
		{
			int up = tmpOrderOut - 3 - i;
			bottomup[up] = tmpOrder[i];
			bottomup[up+1] = tmpOrder[i+1];
			bottomup[up+2] = tmpOrder[i+2];
		}
		bottomup[tmpOrderOut] = -1;
	}
	

	
}
