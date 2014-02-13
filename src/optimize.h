#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include "character.h"

namespace network { struct data; }

namespace optimize
{
	struct optstate;

	optstate *create(int characters, int allocnodes);
	void free(optstate *s);
	void init();
	void copy(optstate *target, optstate *source);
	
	// returns new distance
	character::distance_t optimize(network::data *data, bool all_chars = false);
}

#endif
