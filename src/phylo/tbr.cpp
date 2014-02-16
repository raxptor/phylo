#include "tbr.h"
#include "network.h"
#include "newick.h"
#include "character.h"
#include "optimize.h"
#include <cstring>
#include <iostream>
#include <mtw/mersenne-twister.h>
#include <cstdlib>

namespace tbr
{	
	//#define DPRINT(x) { std::cout << x << std::endl; }
	#define DPRINT(x) { }

	// This will recompute the results	
	//
	// #define CHECK_RESULTS
	//#define CHECK_ALL_BEST

	enum
	{
		MAX_NODES	= 4096
	};
	
	long long networks = 0;
	
	inline void count_networks()
	{
		if (++networks % 1000000 == 0)
			std::cout << "[tbr] - searched " << (networks/1000000) << "M networks.." << std::endl;
	}
	
	void bisect(network::data *d, network::idx_t n0, network::idx_t n1, network::idx_t *s0, network::idx_t *s1)
	{
		// fan outs from n0 and n1 the other way around
		network::idx_t out0[3] = { network::NOT_IN_NETWORK, network::NOT_IN_NETWORK, 0 };
		network::idx_t out1[3] = { network::NOT_IN_NETWORK, network::NOT_IN_NETWORK, 0 };
	
		if (n0 < d->mtx_taxons || n1 < d->mtx_taxons)
		{
			DPRINT("Ignoring single taxon for now.");
			exit(1);
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
			exit(1);
			return;
		}
		
		outs = out1;		
		if (N1->c0 != n0) *outs++ = N1->c0;
		if (N1->c1 != n0) *outs++ = N1->c1;
		if (N1->c2 != n0) *outs++ = N1->c2;
		if (outs != &out1[2])
		{
			DPRINT("Network is corrupted because fanout != 2 from " << n1);
			exit(1);
			return;
		}
		
		network::edge_merge(d->network, out0[0], out0[1], n0);
		network::edge_merge(d->network, out1[0], out1[1], n1);

		DPRINT("Merge freed up " << n0 << " and " << n1);

		network::node_free(d, n0);
		network::node_free(d, n1);
		
		*s0 = out0[0];
		*s1 = out1[0];
	}
	
	void roulette(network::data *d, int org_dist, network::idx_t s0, network::idx_t s1, output *out, bool single_source)
	{ 		
		DPRINT("--- Roulette starts with networks (" << s0 << ") and (" << s1 << ")" << " ---");
		
		network::edgelist net0, net1;
		network::trace_edgelist(d, s0, &net0);

		int length0, length1;
		
		for (int i=0;i<net0.count;i++)
		{
			if (net0.pairs[i] < d->mtx_taxons)
			{
				DPRINT("###### OPTIMIZING FIRST NETWORK " << net0.pairs[i]);
				
				length0 = optimize::optimize(d, net0.pairs[i], true);
		
				// newick::print(d, net0.pairs[i]);
				DPRINT(" length was " << length0);
				break;
			}
		}
				if (!single_source)
		{
			network::trace_edgelist(d, s1, &net1);

			for (int i=0;i<net1.count;i++)
			{
				if (net1.pairs[i] < d->mtx_taxons)
				{
					DPRINT("###### OPTIMIZING SECOND NETWORK " << net1.pairs[i]);
					length1 = optimize::optimize(d, net1.pairs[i], true);
					// newick::print(d, net1.pairs[i]);
					DPRINT("length was " << length1);
					break;
				}
			}
		}
		else
		{
			// is taxon
			length1 = 0;
			net1.count = 2;
			net1.pairs[0] = s1;
			net1.pairs[1] = network::NOT_IN_NETWORK;
			DPRINT("#### SOURCE IS SINGLE (" << s1 << ") LENGTH = 0");
		}		
		
/*
		DPRINT("--- Left nework --- (d=" << length0 << ")");
		for (int i=0;i<net0.count;i++)
		{
			if (net0.pairs[i] < d->mtx_taxons)
			{
				newick::print(d, net0.pairs[i]);
				break;
			}
		}
		
		DPRINT("--- Right nework --- (d=" << length1 << ")");
		for (int i=0;i<net1.count;i++)
		{
			if (net1.pairs[i] < d->mtx_taxons)
			{
				newick::print(d, net1.pairs[i]);
				break;
			}
		}
*/
		
		// diff between original tree and the sum of clipped.
		network::idx_t a = network::node_alloc(d); // store the computation there
		network::idx_t b = !single_source ? network::node_alloc(d) : -1;

		for (int i=0;i<net0.count;i+=2)
		{
			const network::idx_t _a0 = net0.pairs[i];
			const network::idx_t _a1 = net0.pairs[i+1];
			const network::node a0 = d->network[_a0];
			const network::node a1 = d->network[_a1];
			

			// alright, this outer is the source tree now
			//
			DPRINT("Prep source tree " << _a0 << " and " << _a1 << " to " << a);
			optimize::prepare_source_tree_root(d, _a0, _a1, a);


			for (int j=0;j<net1.count;j+=2)
			{
				const network::idx_t _b0 = net1.pairs[j];
				const network::idx_t _b1 = net1.pairs[j+1];
			
				// At this point _b1 can potentially be -1 (NOT_IN_NETWORK)
				DPRINT("Gluing together " << _a0 << "-" << _a1 << " with " << _b0 << "-" << _b1 <<  " tmp nodes " << a << " and " << b);
				
				// 
				const int mergediff = optimize::clip_merge_dist(d, a, _b0, _b1);
				const int newlength = (length0 + length1 + mergediff);
				
				DPRINT("merge diff = " << mergediff << " newlength = " << newlength << " prevlength=" << out->length);
				
				char treebuf[4096];
				
				if (newlength < out->length)
				{
					out->length = newlength;
					network::data *target = out->best_network;
					network::copy(target, d);
					
					if (!single_source)
					{
						network::edge_split(target->network, _a0, _a1, a);
						network::edge_split(target->network, _b0, _b1, b);
						target->network[a].c2 = b;
						target->network[b].c2 = a;
					}
					else
					{
						// we copy the allocated tmp vertex now which will be needed for the insert
						// so free that one back
						network::node_free(target, a);
						// taxon rejoined
						DPRINT("inserting " << _b0 << " to " << _a0 << " " << _a1 << " tot " << target->allocnodes);
						network::insert(target, _a0, _a1, _b0);
					}
					
					#if defined(CHECK_ALL_BEST) || defined(CHECK_RESULTS)
					int real = optimize::optimize(target);
					if (real != newlength)
					{
						std::cout << "TBR found false length when copying! nl=" << newlength << " but " << real << std::endl;
						optimize::print_state(d->opt, d->allocnodes, d->allocnodes);
						exit(0);
					}
					#endif
					
					network::sort(target);
					char *outbuf = treebuf;
					network::to_string_2(target, &outbuf, 0, network::NOT_IN_NETWORK);
					*outbuf = 0;
					
					DPRINT("sorted is " << treebuf);
					
					out->equal_length.clear();
					out->equal_length.insert(treebuf);
				}
				else if (newlength == out->length)
				{
					// tmp construct to make string	
					const network::node b0 = d->network[_b0];
					const network::node b1 = d->network[_b1 >= 0 ? _b1 : 0];

					network::edge_split(d->network, _a0, _a1, a);
					
					if (!single_source)
					{
						network::edge_split(d->network, _b0, _b1, b);
						d->network[a].c2 = b;
						d->network[b].c2 = a;
					}
					else
					{
						d->network[a].c2 = _b0;
					}
					
					char *bufptr = treebuf;
					network::sort(d);
					network::to_string_2(d, &bufptr, 0, network::NOT_IN_NETWORK);
					*bufptr = 0;
					
					DPRINT("sorted is " << treebuf);
					
					out->equal_length.insert(treebuf);
					
					if (rand_u32()%3==0)
						network::copy(out->best_network, d);

					d->network[_b0] = b0;
					
					if (!single_source)
						d->network[_b1] = b1;
					
					d->network[_a0] = a0;
					d->network[_a1] = a1;
				}
				
				#if defined(CHECK_RESULTS)

				if (!single_source)
				{
					// verify	
					network::data *tmp = network::alloc(d->matrix);
					network::copy(tmp, d);
					network::edge_split(tmp->network, _a0, _a1, a);
					network::edge_split(tmp->network, _b0, _b1, b);
					tmp->network[a].c2 = b;
					tmp->network[b].c2 = a;

					const int actual = optimize::optimize(tmp);

					DPRINT(" => This would result in " << mergediff << " extra L0 & L1 =" << length0 << " " << length1);
					DPRINT(" => Computed length = " << (length0 + length1 + mergediff));
					DPRINT(" => Actual length = " << actual);

					if (actual != newlength)
					{
						optimize::print_state(d->opt, d->allocnodes, d->allocnodes);
						std::cerr << "Optimize computed wrong merge length, " << newlength << " vs " << actual << std::endl;
						exit(-2);
					} 
					
					network::free(tmp);
				}
				#endif

				count_networks();
			}
		}
		
		if (a != -1)
			network::node_free(d, a);
		if (b != -1)
			network::node_free(d, b);
	}

	//
	int run(network::data *d, output * out)
	{

		network::edgelist tmp;
		network::trace_edgelist(d, 0, &tmp);
		character::distance_t org_dist = optimize::optimize(d);

		DPRINT("tbr with " << tmp.count/2 << " edges to bisect, distance = " << org_dist);
		
		network::data *bisected = network::alloc(d->matrix);

		for (int i=0;i<tmp.count;i+=2)
		{
			if (tmp.pairs[i] >= d->mtx_taxons && tmp.pairs[i+1] >= d->mtx_taxons)
			{
				// lets work on a temp copy since bisecting and rouletting leaves it bisected.
				// could maybe join them together the original way for the last step
				network::copy(bisected, d);
				network::idx_t s0, s1;

				bisect(bisected, tmp.pairs[i], tmp.pairs[i+1], &s0, &s1);
			
				roulette(bisected, org_dist, s0, s1, out, false);

			}
			else
			{
				network::copy(bisected, d);
				network::idx_t taxon = tmp.pairs[i];
				network::idx_t inner = tmp.pairs[i+1];
				if (taxon >= d->mtx_taxons)
				{
					taxon = tmp.pairs[i+1];
					inner = tmp.pairs[i];
				}
			
				// remember where it goes and insert it.
				network::idx_t r0, r1;
				network::the_two_others(bisected->network, inner, taxon, &r0, &r1);				
				DPRINT("i am (" << taxon << ") others are " << r0 << "," << r1);
				
				network::disconnect(bisected, taxon);
				
				
				roulette(bisected, org_dist, r0, taxon, out, true);
				
				
			}
		}
		
		network::free(bisected);
		return out->length < org_dist;
	}
}
