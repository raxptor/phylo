#include "ratchet.h"
#include "character.h"
#include "matrix.h"
#include "network.h"
#include "newick.h"
#include "tbr.h"

#include <cstdlib>
#include <iostream>

#define DPRINT(x) { std::cout << x << std::endl; }

namespace ratchet
{
	void run(network::data *d)
	{
		int boosts = d->mtx_characters * 3 / 10 + 1;
		
		int picks[1024];
		for (int i=0;i<boosts;i++)
			picks[i] = rand() % d->mtx_characters;

		const int old = d->dist;
		
		network::data *new_net = network::alloc(d->matrix);
		network::copy(new_net, d);
		
		// manipulate & recalc
		for (int i=0;i<d->allocnodes;i++)
			character::toggle_boost(new_net->characters[i], picks, boosts);
		network::recompute_dist(new_net);	
			
		// run
		for (int i=0;i<100;i++)
			if (!tbr::run(new_net))
				break;

		// un-do boost
		for (int i=0;i<d->allocnodes;i++)
			character::toggle_boost(new_net->characters[i], picks, boosts);	
		network::recompute_dist(new_net);	
		
		if (new_net->dist < old)
		{
			std::cout << "ratchet: found net (ph1) with dist " << new_net->dist << " :";
			newick::print(new_net);
		}
		
		// run again
		for (int i=0;i<100;i++)
			if (!tbr::run(new_net))
				break;

		if (new_net->dist < old)
		{
			std::cout << "ratchet: found net (ph2) with dist " << new_net->dist << " :";
			newick::print(new_net);
		}
		
		network::copy(d, new_net);
		network::free(new_net);
	}
}
