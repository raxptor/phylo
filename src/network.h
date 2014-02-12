#ifndef __NETWORK_H__
#define __NETWORK_H__

#include "matrix.h"
#include "character.h"

namespace network
{
	typedef int idx_t;
	
	enum {
		NOT_IN_NETWORK = -1,
		MAX_NODES = 2048
	};
	
	struct node
	{
		idx_t c0, c1, c2;
	};

	struct edgelist
	{
		int count;
		network::idx_t pairs[2*MAX_NODES];
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
		character::distance_t *dist_cache;
	};
	
	data* alloc(matrix::data *mtx);
	void copy(data *target, data *source);
	void recompute_dist(data *target);
	
	void init(data *d, idx_t taxon0, idx_t taxon1);
	
	idx_t insert(data *d, idx_t n0, idx_t n1, idx_t which);
	void disconnect(data *d, idx_t taxon);
	
	idx_t node_alloc(data *d);
	void node_free(data *d, idx_t where);

	// splits edge n0,n1 and inserts newmiddle, pointing its c0 and c1 to n0 and n1
	void edge_split(node *net, idx_t n0, idx_t n1, idx_t newmiddle);
	void edge_merge(node *net, idx_t n0, idx_t n1, idx_t middle);
	void the_two_others(node *network, idx_t where, idx_t which, idx_t *r0, idx_t *r1);
	void trace_edgelist(data *d, idx_t start, edgelist * out);
	
	// make all p0 point to the root 
	void treeify(network::data *data, idx_t root);
	
	void print_characters(data *d);

	bool character_fiddle(data *, int target, int tries);
	
	void free(data *d);
}

#endif
