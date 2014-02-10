#include "tree.h"
#include <cstdlib>
#include <iostream>

namespace tree
{
	data* alloc(matrix::data *mtx)
	{
		data *d = new data();
		d->nodes = 0;
		d->matrix = mtx;
		d->allocnodes = 2 * mtx->taxons - 2;
		d->tree = new node[d->allocnodes];
		d->extra_ptr = d->allocnodes - 1;
		return d;
	}

	void init(tree::data *d, int taxon0, int taxon1)
	{
		int it = 0;
		int nit = 0;
		for (unsigned int i=0;i<d->allocnodes;i++)
		{
			d->tree[i].parent = NOT_IN_TREE;
		}
	
		d->tree[taxon0].parent = NONE;
		d->tree[taxon0].sibling = NONE;
		d->tree[taxon1].parent = taxon0;
		d->tree[taxon1].sibling = NONE;
		d->nodes = 2;
		d->extra_ptr = d->allocnodes -1;
	}
	
	idx_t insert(tree::data *d, int where, int taxon)
	{
		//
		const int p = d->tree[where].parent;
		if (p < 0)
		{
			std::cerr << "err: insert on target node with no parent " << where << std::endl;
			exit(-1);
		}
		
		if (d->tree[taxon].parent != NOT_IN_TREE)
		{
			std::cerr << "err: taxon is already in tree" << std::endl; 
			exit(-1);
		}
		
		idx_t n = alloc(d);
		
		d->tree[where].parent = n;
		d->tree[where].sibling = taxon;
		d->tree[n].parent = p;
		d->tree[n].sibling = NONE;
		d->tree[taxon].parent = n;
		d->tree[taxon].sibling = where;
		return n;
	}
	
	int alloc(data *d)
	{
		idx_t e = d->extra_ptr;
			
		for (int i=0;i<d->allocnodes;i++)
		{
			if (++e == d->allocnodes)
				e = d->matrix->taxons;
		
			if (d->tree[e].parent == NOT_IN_TREE)
			{
				d->extra_ptr = e;
				return e;
			}
		}
		
		std::cerr << "[HELP] extra alloc failed" << std::endl;
		exit(-2);
		
		d->extra_ptr = e;
		return -1;
	}
	
	// this is O(n^2)!
	std::string newick_subtree(data *d, int where)
	{
		std::string cont("(");
		bool first = true;

		for (int i=0;i<d->allocnodes;i++)
		{
			if (d->tree[i].parent == where)
			{
				if (!first) 
					cont.append(",");
				else
					first = false;
					
				cont.append(newick_subtree(d, i));
			}
		}
		
		cont.append(")");
		
		if (where < d->matrix->taxons)
			cont.append(matrix::taxon_name(d->matrix, where));

		return cont;		
	}

	void print_newick(data *d)
	{
		std::cout << "Newick: " << newick_subtree(d, 0) << ";" << std::endl;
	}
	
	void free(data *d)
	{
		delete [] d->tree;
		delete d;
	}
}
