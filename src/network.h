#ifndef __NETWORK_H__
#define __NETWORK_H__

#include "matrix.h"
#include "character.h"

namespace network
{
	typedef short idx_t;
	
	enum {
		NOT_IN_NETWORK = -1
	};
	
	struct node
	{
		idx_t c0, c1, c2;
	};
	
	struct idata;
	struct data
	{
		int allocnodes;
		int mtx_taxons, mtx_characters;

		matrix::data *matrix;
		node *network;
		
		character::state_t **characters;
		character::buf *cbuf;
		
		int dist;
		int freecount;
		idx_t *freelist;
	};
	
	data* alloc(matrix::data *mtx);
	
	void init(data *d, idx_t taxon0, idx_t taxon1);

	idx_t insert(data *d, idx_t n0, idx_t n1, idx_t taxon);
	void disconnect(data *d, idx_t taxon);

	idx_t node_alloc(data *d);
	void node_free(data *d, idx_t where);
	
	void free(data *d);
}

#endif
