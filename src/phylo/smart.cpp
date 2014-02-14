#include "smart.h"
#include "network.h"
#include "matrix.h"
#include "optimize.h"

#include <vector>
#include <cstdlib>
#include <iostream>

#include <mtw/mersenne-twister.h>

namespace smart
{
	network::data *make(matrix::data *mtx)
	{
		if (mtx->taxons < 2)
		{
			std::cerr << "Trying to make a network with less than 2 taxons!" << std::endl;
			return 0;
		}
		
		network::data *t = network::alloc(mtx);
		make_inplace(t);
		return t;
	}
	
	struct edge
	{
		network::idx_t n0, n1;
	};
	
	void make_inplace(network::data *t)
	{
		std::vector<network::idx_t> taxons;
		
		for (int i=0;i<t->matrix->taxons;i++)
			taxons.push_back(i);
		
		const int p1 = rand_u32() % taxons.size();
		const int t1 = taxons[p1];
		taxons.erase(taxons.begin() + p1);
		const int p2 = rand_u32() % taxons.size();
		const int t2 = taxons[p2];
		taxons.erase(taxons.begin() + p2);
		
		network::init(t, t1, t2);
		
		while (!taxons.empty())
		{
			const int pick = rand_u32() % taxons.size();
		
			network::edgelist out;
			network::trace_edgelist(t, t1, &out);
			
			int select = 0;
			int min = 100000000;
			for (int i=0;i<out.count;i+=2)
			{
				network::idx_t neu = network::insert(t, out.pairs[i], out.pairs[i+1], taxons[pick]);
				int d = optimize::optimize(t, t1);
				if (d < min)
				{
					min = d;
					select = i;
				}
				network::disconnect(t, taxons[pick]);
			}
			
			network::insert(t, out.pairs[select], out.pairs[select+1], taxons[pick]);
			
			taxons.erase(taxons.begin() + pick);
		}
	}
}
