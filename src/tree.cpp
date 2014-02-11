#include "tree.h"
#include "character.h"

#include <cstdlib>
#include <iostream>
#include <vector>

//#define TREEDEBUG
//#define TREECHECKS

#if defined(TREEDEBUG)
	#define DPRINT(x) { std::cout << x << std::endl; }
#else
	#define DPRINT(x) { }
#endif

namespace tree
{
	data* alloc(matrix::data *mtx)
	{
		data *d = new data();
		d->matrix = mtx;
		d->mtx_taxons = mtx->taxons;
		d->mtx_characters = mtx->characters;
		d->allocnodes = 2 * mtx->taxons - 2;
		d->tree = new node[d->allocnodes];
		d->characters = new character::state_t*[d->allocnodes];
		d->cbuf = character::alloc(mtx->characters, d->allocnodes, d->characters);
		
		for (int i=0;i<mtx->taxons;i++)
			character::copy(d->characters[i], mtx->taxonbase[i], mtx->characters);

		// free list
		d->dist = 0;
		d->root = 0;
		d->freecount = d->allocnodes - mtx->taxons;
		d->freelist = new idx_t[d->freecount];
		for (unsigned int i=0;i<d->freecount;i++)
			d->freelist[i] = i + mtx->taxons;
			
		return d;
	}

	void free(data *d)
	{
		delete [] d->characters;
		delete [] d->tree;
		character::free(d->cbuf);
		delete d;
	}

	void init(tree::data *d, idx_t taxon0, idx_t taxon1)
	{
		for (unsigned int i=0;i<d->allocnodes;i++)
			d->tree[i].parent = NOT_IN_TREE;

		d->freecount = d->allocnodes - d->matrix->taxons;
		for (unsigned int i=0;i<d->freecount;i++)
			d->freelist[i] = i + d->matrix->taxons;
	
		d->dist = character::distance(d->characters[taxon0], d->characters[taxon1], d->mtx_characters);
		d->tree[taxon0].parent = NONE;
		d->tree[taxon0].sibling = NONE;
		d->tree[taxon1].parent = taxon0;
		d->tree[taxon1].sibling = NONE;
		d->root = taxon0;
		
		DPRINT("Tree starts with " << taxon0 << " and " << taxon1 << std::endl);
	}

	void check(tree::data *d)
	{
#if defined(TREEDEBUG)
		int sum = 0;
		int nodes = 0;	
		for (int i=0;i<d->allocnodes;i++)
		{
			idx_t sib = d->tree[i].sibling;
			idx_t parent = d->tree[i].parent;
			
			if (parent == NOT_IN_TREE)
				continue;
				
			nodes++;
			
			if (parent >= 0 && d->tree[parent].parent != NOT_IN_TREE)
				sum += character::distance(d->characters[i], d->characters[parent], d->mtx_characters);
			
			if (sib >= 0)
			{
				if (sib >= d->allocnodes)
				{
					std::cout << i << " has sibling " << sib << " out of range!" << std::endl;
					continue;
				}
					
				if (d->tree[sib].parent == NOT_IN_TREE)
				{
					std::cout << i << " has sibling (" << sib << " which does not exist" << std::endl;
				}
				else if (d->tree[sib].parent != NONE)
				{
					if (d->tree[sib].sibling != i)
						std::cout << i << " has sibling " << sib << " which points to " << d->tree[sib].sibling << std::endl;
					if (d->tree[sib].parent != parent)
						std::cout << i << " has sibling " << sib << " with a different parent " << d->tree[sib].parent << " should be " << parent << std::endl;
				}
			}
			
			if (d->tree[parent].parent == NOT_IN_TREE)
			{
				std::cout << i << " has parent " << parent << " which does not exist in the tree " << std::endl;
			}
		}
		
		if (sum != d->dist)
			std::cout << "Tree distance computed to " << sum << " but tree says " << d->dist << " nodes=" << nodes << std::endl;
#endif
	}
	
	idx_t insert(tree::data *d, idx_t where, idx_t taxon)
	{
		const int par = d->tree[where].parent;
		const idx_t sib = d->tree[where].sibling;

#if defined(TREECHECKS)
		if (par < 0)
		{
			std::cerr << "err: insert " << taxon << " on target node " << where << " with no parent " << where << std::endl;
			exit(-1);
		}
		
		if (d->tree[taxon].parent != NOT_IN_TREE)
		{
			std::cerr << "err: taxon (" << taxon << ") is already in tree, parent=" << d->tree[taxon].parent << std::endl; 
			exit(-1);
		}
#endif
		
		const idx_t n = node_alloc(d);
		const int chars = d->mtx_characters;
		

		DPRINT(" insert(" << taxon << ") @ " << where << " newnode=" << n);
		
		if (sib >= 0)
		{
			d->tree[sib].sibling = n;
		}
		
		d->tree[where].parent = n;
		d->tree[where].sibling = taxon;
		d->tree[n].parent = par;
		d->tree[n].sibling = sib;
		d->tree[taxon].parent = n;
		d->tree[taxon].sibling = where;

		// generate 
		character::threesome(d->characters[par], d->characters[where], d->characters[taxon], d->characters[n], d->mtx_characters);

		// maybe all of this could be merged into one calculation with the above
		d->dist -= character::distance(d->characters[where], d->characters[par], chars);
		d->dist += character::distance(d->characters[where], d->characters[n], chars);
		d->dist += character::distance(d->characters[n], d->characters[par], chars);
		d->dist += character::distance(d->characters[n], d->characters[taxon], chars);

		DPRINT("     => d = " << d->dist);

#if defined(TREECHECKS)		
		if (sib >= 0 && d->tree[sib].parent != d->tree[n].parent)
		{
			std::cerr << " insert bugged out " << std::endl;
		}
		check(d);
#endif
		return n;
	}
	
	void disconnect(tree::data *d, idx_t which)
	{
#if defined(TREECHECKS)
		if (d->tree[which].parent < 0)
		{
			std::cerr << "err: cannot erase node " << which << std::endl;
			return;
		}
		
		if (d->tree[d->tree[which].parent].parent == NOT_IN_TREE)
		{
			std::cerr << "err: disconnecting node with parent not in tree (" << d->tree[which].parent << ")" << std::endl;
			return;
		}
#endif
		
		const idx_t sib = d->tree[which].sibling;
		const idx_t par = d->tree[which].parent;
		const idx_t parsib = d->tree[par].sibling;
		const idx_t parpar = d->tree[par].parent;
	
		// uncomment for some tree debugging	
		DPRINT("Erasing node " << which << " par:" << par << " sib:" << sib << " parsib:" << parsib << " parpar:" << parpar);
		
#if defined(TREECHECKS)
		if (sib < 0)
			std::cout << "no sibling here?" << std::endl;
#endif

		d->tree[sib].parent = parpar;
		d->tree[sib].sibling = parsib;
		
		if (parsib >= 0)
		{
			d->tree[parsib].sibling = sib;
		}

		// which must a terminal node
		d->tree[par].parent = NOT_IN_TREE;
		d->tree[which].parent = NOT_IN_TREE;
		node_free(d, par);
		
		const int chars = d->mtx_characters;
		d->dist -= character::distance(d->characters[which], d->characters[par], chars);
		d->dist -= character::distance(d->characters[par], d->characters[sib], chars);
		d->dist -= character::distance(d->characters[par], d->characters[parpar], chars);
		d->dist += character::distance(d->characters[sib], d->characters[parpar], chars);
		DPRINT("   d => " << d->dist);	

#if defined(TREECHECKS)
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
