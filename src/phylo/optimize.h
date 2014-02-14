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
	
	// returns new distance
	character::distance_t optimize(network::data *data, int root=0);
}

#endif
