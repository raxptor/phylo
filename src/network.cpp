#include "network.h"
#include "character.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cstring>

// Debug stuffs
// #define NETWORKDEBUG
// #define NETWORKCHECKS

#if defined(NETWORKDEBUG)
	#define DPRINT(x) { std::cout << x << std::endl; }
#else
	#define DPRINT(x) { }
#endif

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
		d->characters = new character::state_t*[d->allocnodes];
		d->cbuf = character::alloc(mtx->characters, d->allocnodes, d->characters);
		
		for (int i=0;i<mtx->taxons;i++)
			character::copy(d->characters[i], mtx->taxonbase[i], mtx->characters);

		// free list
		d->dist = 0;
		d->freecount = d->allocnodes - mtx->taxons;
		d->freelist = new idx_t[d->freecount];
		for (unsigned int i=0;i<d->freecount;i++)
			d->freelist[i] = i + mtx->taxons;
			
		return d;
	}
	
	// 
	void copy(data *target, data *source)
	{
		memcpy(target->network, source->network, sizeof(node) * source->allocnodes);
		character::copy(target->cbuf, source->cbuf);
		memcpy(target->freelist, source->freelist, sizeof(idx_t) * target->freecount);
		target->freecount = source->freecount;
		target->dist = source->dist;
	}
	
	void recompute_dist(data *target)
	{
		edgelist e;
		trace_edgelist(target, 0, &e);
		target->dist = 0;
		for (int i=0;i<e.count;i+=2)
		{
			character::distance_t dist = character::distance(target->characters[e.pairs[i]], target->characters[e.pairs[i+1]], target->mtx_characters);
			target->dist += dist;
			DPRINT(" -> recompute dist, dist=" << dist << " for " << e.pairs[i] << "-" << e.pairs[i+1]);
		}
	}

	void free(data *d)
	{
		delete [] d->characters;
		delete [] d->network;
		character::free(d->cbuf);
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
	
		d->dist = character::distance(d->characters[taxon0], d->characters[taxon1], d->mtx_characters);
		d->network[taxon0].c0 = taxon1;
		d->network[taxon1].c0 = taxon0;
		
		DPRINT("Network starts with " << taxon0 << " and " << taxon1 << std::endl);
	}

	void check(network::data *d)
	{
#if defined(NETWORKDEBUG)
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
		const int chars = d->mtx_characters;

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
		
		// generate 
		character::threesome(d->characters[n0], d->characters[n1], d->characters[which], d->characters[n], d->mtx_characters);

		// maybe all of this could be merged into one calculation with the above
		d->dist -= character::distance(d->characters[n0], d->characters[n1], chars);
		d->dist += character::distance(d->characters[n0], d->characters[n], chars);
		d->dist += character::distance(d->characters[n], d->characters[n1], chars);
		d->dist += character::distance(d->characters[n], d->characters[which], chars);

		DPRINT("     => d = " << d->dist << " new=" << n);

#if defined(NETWORKCHECKS)		
		check(d);
#endif

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
				
		//
		idx_t r0, r1;
		the_two_others(network, n, which, &r0, &r1);
		
		DPRINT(" -> merging from " << r0 << "-" << r1 << " with " << n << " in the middle");
		edge_merge(d->network, r0, r1, n);		
			
		// which must a terminal node
		node_free(d, n);
		
		const int chars = d->mtx_characters;
		d->dist -= character::distance(d->characters[which], d->characters[n], chars);
		d->dist -= character::distance(d->characters[r0], d->characters[n], chars);
		d->dist -= character::distance(d->characters[r1], d->characters[n], chars);
		d->dist += character::distance(d->characters[r0], d->characters[r1], chars);
		DPRINT("   d => " << d->dist);	

#if defined(NETWORKCHECKS)
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

	
}
