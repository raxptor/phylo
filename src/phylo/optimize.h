#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include "character.h"

namespace network { struct data; struct node; }
namespace matrix { struct data; }

namespace optimize
{
	enum
	{
		BUFSIZE  = 65768, // char group * taxa limit
		MAX_CHARACTERS = 2048,
		MAX_NODES = 2048
	};
	
	// do not touch these unless you rewrite all SSE code
	#define BLOCKSHIFT 4
	#define BLOCKSIZE (1<<BLOCKSHIFT)
	
	typedef unsigned char st_t;
	
	struct cgroup_data
	{
		int count;
		int taxonwidth;
		int bufsize;
		
		st_t ostate[BUFSIZE] __attribute__ ((aligned(BLOCKSIZE)));
		st_t pstate[BUFSIZE] __attribute__ ((aligned(BLOCKSIZE)));
		st_t fstate[BUFSIZE] __attribute__ ((aligned(BLOCKSIZE)));
		st_t weights[MAX_CHARACTERS] __attribute__ ((aligned(BLOCKSIZE)));
	};
	
	struct optstate
	{
		int root, maxnodes;
		int *first_pass_order;
		
		cgroup_data unordered;
		cgroup_data ordered;
		
		network::node *net;
		matrix::data *matrix;
	};

	optstate *create(network::data *, bool will_copy);
	void free(optstate *s);
	void init();
	void copy(optstate *target, optstate *source);
	void print_state(optstate *st, int maxnodes, int taxons);
	
	// returns new distance
	character::distance_t optimize(network::data *data, int root=0, bool write_final=false);
	
	void reset_taxons_states(cgroup_data *cd, int *first_pass_order, int taxons, int root, int root_htu);
	
	// TBR utilitiy functions
	void prepare_source_tree_root(network::data *d, int s0, int s1, int new_node);
	character::distance_t clip_merge_dist(network::data *d, int source_root, int t0, int t1, int max);

	// TBR ronquist
	void tbr_source_reoptimization(network::data *bisected, network::data *source, int from_where);
 	void tbr_target_reoptimization(network::data *data, network::data *ref, int calcroot, int start);
      
	void set_weight(optstate *st, int pos, int weight);
	void ultranode(network::data *d, int node);
}

#endif
