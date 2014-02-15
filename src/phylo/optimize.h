#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include "character.h"

namespace network { struct data; }

namespace optimize
{
	struct optstate;

	optstate *create(network::data *);
	void free(optstate *s);
	void init();
	void copy(optstate *target, optstate *source);
	void print_state(optstate *st, int maxnodes, int taxons);
	
	// returns new distance
	character::distance_t optimize(network::data *data, int root=0, bool write_final=false);
	character::distance_t clip_merge_dist(network::data *d, int t0, int t1, int s0, int s1);

	// TBR utilitiy functions
	void prepare_source_tree_root(network::data *d, int s0, int s1, int new_node);
	character::distance_t clip_merge_dist(network::data *d, int source_root, int t0, int t1);
}

#endif
