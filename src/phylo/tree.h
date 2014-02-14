#ifndef __TREE_H__
#define __TREE_H__

#include "matrix.h"
#include "character.h"

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
		character::distance_t sum;
	};
	
	struct idata;
	struct data
	{
		int allocnodes;
		int mtx_taxons, mtx_characters;

		matrix::data *matrix;
		node *tree;
		
		character::state_t **characters;
		character::buf *cbuf;
		
		int dist;
		int freecount;
		idx_t root;
		idx_t *freelist;
	};
	
	data* alloc(matrix::data *mtx);
	
	void init(data *d, idx_t taxon0, idx_t taxon1);

	idx_t insert(data *d, idx_t where, idx_t taxon);
	void disconnect(data *d, idx_t which);

	idx_t node_alloc(data *d);
	void node_free(data *d, idx_t where);
	
	void free(data *d);
}

#endif
