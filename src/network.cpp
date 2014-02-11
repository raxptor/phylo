#include "network.h"
#include "character.h"

#include <cstdlib>
#include <iostream>
#include <vector>

//#define NETWORKDEBUG
//#define NETWORKCHECKS

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
	
	idx_t insert(data *d, idx_t n0, idx_t n1, idx_t taxon)
	{
		DPRINT(" Insert " << taxon << " at " << n0 << "-" << n1);
		
		node * const net = d->network;
		
		const idx_t n = node_alloc(d);
		const int chars = d->mtx_characters;

		// insert middle		
		edge_split(d->network, n0, n1, n);

		net[n].c2 = taxon;		
		net[taxon].c0 = n;
		net[taxon].c1 = NOT_IN_NETWORK;
		net[taxon].c2 = NOT_IN_NETWORK;
		
		// generate 
		character::threesome(d->characters[n0], d->characters[n1], d->characters[taxon], d->characters[n], d->mtx_characters);

		// maybe all of this could be merged into one calculation with the above
		d->dist -= character::distance(d->characters[n0], d->characters[n1], chars);
		d->dist += character::distance(d->characters[n0], d->characters[n], chars);
		d->dist += character::distance(d->characters[n], d->characters[n1], chars);
		d->dist += character::distance(d->characters[n], d->characters[taxon], chars);

		DPRINT("     => d = " << d->dist << " new=" << n);

#if defined(NETWORKCHECKS)		
		check(d);
#endif

		return n;
	}
	
	void disconnect(data *d, idx_t taxon)
	{
		node * network = d->network;

		const idx_t n = network[taxon].c0;
		
		//
		idx_t r0, r1;
		if (network[n].c0 == taxon)
		{
			r0 = network[n].c1;
			r1 = network[n].c2;
		}
		else if (network[n].c1 == taxon)
		{
			r0 = network[n].c0;
			r1 = network[n].c2;
		}
		else if (network[n].c2 == taxon)
		{
			r0 = network[n].c0;
			r1 = network[n].c1;
		}
		
		edge_merge(d->network, r0, r1, n);		
			
		// which must a terminal node
		node_free(d, n);
		
		const int chars = d->mtx_characters;
		d->dist -= character::distance(d->characters[taxon], d->characters[n], chars);
		d->dist -= character::distance(d->characters[r0], d->characters[n], chars);
		d->dist -= character::distance(d->characters[r1], d->characters[n], chars);
		d->dist += character::distance(d->characters[r0], d->characters[r1], chars);
		DPRINT("   d => " << d->dist);	

#if defined(NETWORKCHECKS)
		check(d);
#endif
	}
	
	inline idx_t node_alloc(data *d)
	{
		if (!d->freecount)
		{
			std::cerr << "[HELP] extra alloc failed" << std::endl;
			exit(-2);
		}
		
		return d->freelist[--d->freecount];
	}
	
	inline void node_free(data *d, idx_t where)
	{
		if (where < d->matrix->taxons)
		{
			std::cerr << "Trying to free original taxon" << std::endl;
		}
		
		d->freelist[d->freecount++] = where;
	}
	
}
