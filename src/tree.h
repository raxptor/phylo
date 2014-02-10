#ifndef __TREE_H__
#define __TREE_H__

#include "matrix.h"

namespace tree
{
	typedef int idx_t;
	
	enum {
		NONE = -1,
		NOT_IN_TREE = -2
	};
	
	struct node
	{
		idx_t parent;
		idx_t sibling;
	};
	
	struct data
	{
		int nodes, allocnodes;
		matrix::data *matrix;
		node *tree;
		idx_t *in_tree;
		idx_t *not_in_tree;
		idx_t extra_ptr;
	};
	
	data* alloc(matrix::data *mtx);
	
	void init(data *d, int taxon0, int taxon1);
	idx_t insert(data *d, int where, int taxon);
	int alloc(data *d);
	
	void print_newick(data *d);
	
	void free(data *d);
}

#endif
