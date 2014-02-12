#include "ratchet.h"
#include "character.h"
#include "matrix.h"
#include "network.h"
#include "newick.h"
#include "tbr.h"

#include <cstdlib>
#include <iostream>
#include <mtw/mersenne-twister.h>

#define DPRINT(x) { std::cout << x << std::endl; }

namespace ratchet
{
	void run(network::data *d, tbr::output *out)
	{
		int boosts = rand() % (d->mtx_characters / 5 + 1) + 1;
		
		int picks[1024];
		for (int i=0;i<boosts;i++)
			picks[i] = rand_u32() % d->mtx_characters;

		const int old = d->dist;
		
		network::data *new_net = network::alloc(d->matrix);
		network::copy(new_net, d);
		
		// manipulate & recalc
		for (int i=0;i<d->allocnodes;i++)
			character::toggle_boost(new_net->characters[i], picks, boosts);
		network::recompute_dist(new_net);	
			
		// set up temporary output structure for tbr, we won't record any of these
		// networks as final		

		tbr::output tout;
		tout.best_network = network::alloc(d->matrix);
		network::copy(tout.best_network, new_net);

		for (int i=0;i<100;i++)
			if (!tbr::run(new_net, &tout))
				break;

		network::free(tout.best_network);

		// un-do boost
		for (int i=0;i<d->allocnodes;i++)
			character::toggle_boost(new_net->characters[i], picks, boosts);
		network::recompute_dist(new_net);	
		
		if (new_net->dist < out->best_network->dist)
		{
			std::cout << "ratchet: found net (ph1) with dist " << new_net->dist << " :";
			network::copy(out->best_network, new_net);
			newick::print(new_net);
		}
		
		const character::distance_t prestore = out->best_network->dist;
		
		// run again (no temp out this time)
		for (int i=0;i<100;i++)
			if (!tbr::run(new_net, out))
				break;

		if (new_net->dist < prestore)
		{
			std::cout << "ratchet: found net (ph2) with dist " << new_net->dist << " :";
			newick::print(new_net);
		}

		// keep worknig
		if (new_net->dist <= prestore)
		{
			network::copy(out->best_network, new_net);
		}
		
		network::copy(d, new_net);
		network::free(new_net);
	}
}
