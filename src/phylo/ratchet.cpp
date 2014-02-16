#include "ratchet.h"
#include "character.h"
#include "matrix.h"
#include "network.h"
#include "newick.h"
#include "tbr.h"
#include "optimize.h"

#include <cstdlib>
#include <iostream>
#include <mtw/mersenne-twister.h>

#define DPRINT(x) { std::cout << x << std::endl; }

namespace ratchet
{
	void run(network::data *d, tbr::output *out)
	{
		int boosts = (d->mtx_characters / 5) + 1;
		int picks[1024];
		for (int i=0;i<boosts;i++)
			picks[i] = rand_u32() % d->mtx_characters;

		const int old = out->length;
		
		network::data *new_net = network::alloc(d->matrix);
		network::copy(new_net, d);

		int nw = 10;
		if (rand_u32()%10 > 5)
			nw = 0;
		else
			nw = 10;
			
		// manipulate
		for (int i=0;i<boosts;i++)
		{
			optimize::set_weight(new_net->opt, picks[i], nw);
		}
		
		// set up temporary output structure for tbr, we won't record any of these
		// networks as final		

		tbr::output tout;
		tout.best_network = network::alloc(d->matrix);
		tout.length = 1000000;
		network::copy(tout.best_network, new_net);

		for (int i=0;i<20;i++)
		{
			if (!tbr::run(new_net, &tout))
				break;
		}

		network::copy(new_net, tout.best_network);
		network::free(tout.best_network);
		
		// restore
		for (int i=0;i<boosts;i++)
			optimize::set_weight(new_net->opt, picks[i], 1);
		
		// see if we did any better than previous best
		const character::distance_t prestore = optimize::optimize(out->best_network);

		// run again (no temp out this time)
		for (int i=0;i<20;i++)
		{
			if (!tbr::run(new_net, out))
				break;
		}

		network::copy(new_net, out->best_network);
		
		if (out->length < old)
		{
			int d = optimize::optimize(out->best_network);
			
			std::cout << "ratchet: found net (ph2) with dist " << d << " recorded(" << out->length << ") previous(" << old << ")" << std::endl;
			newick::print(new_net);
		}

		network::copy(d, new_net);
		network::free(new_net);
	}
}
