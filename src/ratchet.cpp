#include "ratchet.h"
#include "character.h"
#include "matrix.h"
#include "network.h"
#include "newick.h"
#include "tbr.h"
#include "debug.h"
#include "optimize.h"

#include <cstdlib>
#include <iostream>
#include <mtw/mersenne-twister.h>

#define DPRINT(x) { std::cout << x << std::endl; }

namespace ratchet
{
	void run(network::data *d, tbr::output *out)
	{
		int boosts = rand_u32() % (d->mtx_characters / 10 + 2) + 1;
		int picks[1024];
		for (int i=0;i<boosts;i++)
			picks[i] = rand_u32() % d->mtx_characters;

		const int old = optimize::optimize(d);
		
		network::data *new_net = network::alloc(d->matrix);
		network::copy(new_net, d);
		
		// manipulate & recalc
		for (int i=0;i<d->allocnodes;i++)
			character::toggle_boost(new_net->characters[i], picks, boosts);

		// set up temporary output structure for tbr, we won't record any of these
		// networks as final		

		tbr::output tout;
		tout.best_network = network::alloc(d->matrix);
		network::copy(tout.best_network, new_net);

		for (int i=0;i<10;i++)
		{
			if (!tbr::run(new_net, &tout))
				break;
		}

		network::free(tout.best_network);

		// un-do boost
		for (int i=0;i<d->allocnodes;i++)
			character::toggle_boost(new_net->characters[i], picks, boosts);
			
		network::recompute_dist(new_net);	

		// see if we did any better than previous best
		const character::distance_t prestore = optimize::optimize(out->best_network);

		// run again (no temp out this time)
		for (int i=0;i<10;i++)
		{
			if (!tbr::run(new_net, out))
				break;
				
			network::recompute_dist(new_net);	
			if (new_net->dist < prestore)
			{
				int d = new_net->dist;
				new_net->dist = network::distance_by_edges(new_net);;
				std::cout << "ratchet: found net (ph2) with dist " << d << " actual(" << new_net->dist << ") previous(" << prestore << ")" << std::endl;
				newick::print(new_net);
				network::copy(out->best_network, new_net);
				out->equal_length.clear();
			}
		}

		network::copy(d, new_net);
		network::free(new_net);
	}
}
