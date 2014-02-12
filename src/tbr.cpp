#include "tbr.h"
#include "network.h"
#include "newick.h"
#include "character.h"

#include <cstring>
#include <iostream>
#include <mtw/mersenne-twister.h>
#include <cstdlib>

namespace tbr
{	
	//#define DPRINT(x) { std::cout << x << std::endl; }
	#define DPRINT(x) { }

	enum
	{
		MAX_NODES	= 4096
	};
	
	long long networks = 0;
	
	inline void count_networks()
	{
		if (++networks % 10000000 == 0)
			std::cout << "[tbr] - searched " << (networks/1000000) << "M networks.." << std::endl;
	}

	// dist overrides d->dist
	void announce_and_copy(output *out, network::data *best, int new_dist)
	{
		network::copy(out->best_network, best);
		out->best_network->dist = new_dist;
		network::recompute_dist(out->best_network);
		if (out->best_network->dist != new_dist)
		{
			std::cerr << "err: tbr produced false distance" << std::endl;
			exit(-2);
		}
	}
	
	void bisect(network::data *d, network::idx_t n0, network::idx_t n1, network::idx_t *s0, network::idx_t *s1)
	{
		// fan outs from n0 and n1 the other way around
		network::idx_t out0[3] = { network::NOT_IN_NETWORK, network::NOT_IN_NETWORK, 0 };
		network::idx_t out1[3] = { network::NOT_IN_NETWORK, network::NOT_IN_NETWORK, 0 };
	
		if (n0 < d->mtx_taxons || n1 < d->mtx_taxons)
		{
			DPRINT("Ignoring single taxon for now.");
			return;
		}
		
		network::node *N0 = &d->network[n0];
		network::node *N1 = &d->network[n1];
	
		network::idx_t *outs = out0;
		if (N0->c0 != n1) *outs++ = N0->c0;
		if (N0->c1 != n1) *outs++ = N0->c1;
		if (N0->c2 != n1) *outs++ = N0->c2;
		if (outs != &out0[2])
		{
			DPRINT("Network is corrupted because fanout != 2 from " << n0);
			return;
		}
		
		outs = out1;		
		if (N1->c0 != n0) *outs++ = N1->c0;
		if (N1->c1 != n0) *outs++ = N1->c1;
		if (N1->c2 != n0) *outs++ = N1->c2;
		if (outs != &out1[2])
		{
			DPRINT("Network is corrupted because fanout != 2 from " << n1);
			return;
		}
		
		//
		DPRINT("Disconnecting...");
		network::edge_merge(d->network, out0[0], out0[1], n0);
		network::edge_merge(d->network, out1[0], out1[1], n1);
		
		const int chars = d->mtx_characters;

		// subtract the > and < edges

		d->dist -= character::distance(d->characters[n0], d->characters[out0[0]], chars);
		d->dist -= character::distance(d->characters[n0], d->characters[out0[1]], chars);
		d->dist -= character::distance(d->characters[n1], d->characters[out1[0]], chars);
		d->dist -= character::distance(d->characters[n1], d->characters[out1[1]], chars);
		
		// ..and this edge
		d->dist -= character::distance(d->characters[n0], d->characters[n1], chars);
		
		d->dist += character::distance(d->characters[out0[0]], d->characters[out0[1]], chars);
		d->dist += character::distance(d->characters[out1[0]], d->characters[out1[1]], chars);
		
		network::node_free(d, n0);
		network::node_free(d, n1);
		
		*s0 = out0[0];
		*s1 = out1[0];
	}
	
	void roulette(network::data *d, network::idx_t s0, network::idx_t s1, output *out)
	{ 		
		DPRINT("--- Roulette starts with networks (" << s0 << ") and (" << s1 << ")" << " ---");
		network::edgelist net0, net1;
		network::trace_edgelist(d, s0, &net0);

		/*
		DPRINT("--- Left nework ----");
		for (int i=0;i<net0.count;i++)
		{
			if (net0.pairs[i] < d->mtx_taxons)
			{
				newick::print(d, net0.pairs[i]);
				break;
			}
		}
		*/
		network::trace_edgelist(d, s1, &net1);
		
		/*
		DPRINT("--- Right nework ----");
		for (int i=0;i<net1.count;i++)
		{
			if (net1.pairs[i] < d->mtx_taxons)
			{
				newick::print(d, net1.pairs[i]);
				break;
			}
		}
		*/
		const int chars = d->mtx_characters;
		character::state_t **ch = d->characters;

		network::idx_t a = network::node_alloc(d);
		network::idx_t b = network::node_alloc(d);
		
		for (int i=0;i<net0.count;i+=2)
		{
			const network::idx_t _a0 = net0.pairs[i];
			const network::idx_t _a1 = net0.pairs[i+1];
			const network::node a0 = d->network[_a0];
			const network::node a1 = d->network[_a1];
			
			// edge
			const character::distance_t diff = character::distance(ch[_a0], ch[_a1], chars);
			d->dist -= diff;
			
			for (int j=0;j<net1.count;j+=2)
			{
				const network::idx_t _b0 = net1.pairs[j];
				const network::idx_t _b1 = net1.pairs[j+1];
			
				DPRINT("Gluing together " << _a0 << "-" << _a1 << " with " << _b0 << "-" << _b1);

				// calculate what distance would be here
				character::distance_t new_dist = d->dist - character::distance(ch[_b0], ch[_b1], chars);
				
				// -----------------------------------------------------------
				// TODO: Calculating this four-way thing should be easy to do
				//       much more efficient than this.
				
				// outer a0-a1-b0 => a
				character::threesome(ch[_a0], ch[_a1], ch[_b0], ch[a], chars);
				// inner b0-b1-a => b
				character::threesome(ch[_b0], ch[_b1], ch[a], ch[b], chars);
				
				new_dist += character::distance(ch[_a0], ch[a], chars);
				new_dist += character::distance(ch[_a1], ch[a], chars);
				new_dist += character::distance(ch[_b0], ch[b], chars);
				new_dist += character::distance(ch[_b1], ch[b], chars);
				new_dist += character::distance(ch[a], ch[b], chars);
				
				// actually construct it
				DPRINT("New distance " << new_dist << " vs " << out->best_network->dist);
				if (new_dist < out->best_network->dist || (new_dist == out->best_network->dist && rand_u32()%32 == 0))
				{
					const network::node b0 = d->network[_b0];
					const network::node b1 = d->network[_b1];
					
					network::edge_split(d->network, _a0, _a1, a);
					network::edge_split(d->network, _b0, _b1, b);
					
					// merge 
					d->network[a].c2 = b;
					d->network[b].c2 = a;
					
					announce_and_copy(out, d, new_dist);
					
					// restore
					d->network[_b0] = b0;
					d->network[_b1] = b1;
					d->network[_a0] = a0;
					d->network[_a1] = a1;
				}
				
				/*
				for (int i=0;i<d->allocnodes;i++)
				{
					DPRINT("node[" << i << "] " << d->network[i].c0 << "," << d->network[i].c1 << "," << d->network[i].c2);
				}
				*/
				
				count_networks();
			}
			
			// first net's edge put back in, connections restored
			d->dist += diff;
		}
		
		network::node_free(d, a);
		network::node_free(d, b);		
	}

	//
	int run(network::data *d, output * out)
	{
		network::edgelist tmp;
		network::trace_edgelist(d, 0, &tmp);
		
		character::distance_t org_dist = d->dist;
		
		network::data *td = network::alloc(d->matrix);

		DPRINT("tbr with " << tmp.count/2 << " edges to bisect, distance = " << d->dist << " " << out->best_network->dist);
		
		for (int i=0;i<tmp.count;i+=2)
		{
			if (tmp.pairs[i] >= d->mtx_taxons && tmp.pairs[i+1] >= d->mtx_taxons)
			{
				// lets work on a temp copy since bisecting and rouletting leaves it bisected.
				// could maybe join them together the original way for the last step
				network::copy(td, d);
				network::idx_t s0, s1;
				
				bisect(td, tmp.pairs[i], tmp.pairs[i+1], &s0, &s1);
				roulette(td, s0, s1, out);
			}
			else
			{
				network::idx_t taxon = tmp.pairs[i];
				network::idx_t inner = tmp.pairs[i+1];
				if (taxon >= d->mtx_taxons)
				{
					taxon = tmp.pairs[i+1];
					inner = tmp.pairs[i];
				}
			
				// remember where it goes and insert it.				
				network::idx_t r0, r1;
				network::the_two_others(d->network, inner, taxon, &r0, &r1);				
				network::disconnect(d, taxon);

				// cleanup to not find any junk edges
				network::idx_t upper = d->network[taxon].c0;
				d->network[upper].c0 = d->network[upper].c1 = d->network[upper].c2 = network::NOT_IN_NETWORK;
				d->network[taxon] = d->network[upper];
				

				//			
				network::edgelist np;
				network::trace_edgelist(d, (taxon+1) % d->mtx_taxons, &np);
				for (int j=0;j<np.count;j+=2)
				{
					// this is junk from having not filled in the last nodes

					DPRINT("Re-inserting taxon " << np.pairs[j] << " " << np.pairs[j+1]);
					network::insert(d, np.pairs[j], np.pairs[j+1], taxon);
					
					if (d->dist < out->best_network->dist)
						announce_and_copy(out, d, d->dist);
					
					network::disconnect(d, taxon);
					count_networks();
				}
				
				network::insert(d, r0, r1, taxon);
			}
		}
		
		network::copy(d, out->best_network);
		network::free(td);
		
		return out->best_network->dist < org_dist;
	}
}


