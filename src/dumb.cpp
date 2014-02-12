#include "dumb.h"
#include "network.h"
#include "matrix.h"

#include <vector>
#include <cstdlib>
#include <iostream>

#include <mtw/mersenne-twister.h>

namespace dumb
{
	network::data *make(matrix::data *mtx)
	{
		if (mtx->characters < 2)
		{
			std::cerr << "Trying to make a network with less than 2 characters!" << std::endl;
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
		network::init(t, 0, 1);
		std::vector<edge> edges;
		
		edges.reserve(t->matrix->taxons * 2);
		
		edge e;
		e.n0 = 0;
		e.n1 = 1;
		edges.push_back(e); 
		
		for (int i=2;i<t->matrix->taxons;i++)
		{
			const int pick = rand_u32() % edges.size();
			
			edge tmp = edges[pick];
			
			network::idx_t neu = network::insert(t, edges[pick].n0, edges[pick].n1, i);
			
			edges[pick].n0 = tmp.n0;
			edges[pick].n1 = neu;
			
			edge e1;
			e1.n0 = i;
			e1.n1 = neu;
			
			edge e2;
			e2.n0 = neu;
			e2.n1 = tmp.n1;
			
			edges.push_back(e1);
			edges.push_back(e2);
		}
	}
}
