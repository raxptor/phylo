#ifndef __NETWORK_H__
#define __NETWORK_H__

#include "matrix.h"
#include "character.h"

namespace network
{
	typedef int idx_t;
	
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
		node *network;
		character::state_t **characters;
		int mtx_taxons, mtx_characters;
		int dist;
		int freecount;
		idx_t *freelist;

		// never touched during processing
		int allocnodes;
		matrix::data *matrix;		
		character::buf *cbuf;
	};
	
	data* alloc(matrix::data *mtx);
	
	void init(data *d, idx_t taxon0, idx_t taxon1);
	idx_t insert(data *d, idx_t n0, idx_t n1, idx_t which);
	void disconnect(data *d, idx_t taxon);
	
	idx_t node_alloc(data *d);
	void node_free(data *d, idx_t where);

	// splits edge n0,n1 and inserts newmiddle, pointing its c0 and c1 to n0 and n1
	void edge_split(node *net, idx_t n0, idx_t n1, idx_t newmiddle);
	void edge_merge(node *net, idx_t n0, idx_t n1, idx_t middle);
	
	void free(data *d);
}

#endif
