
#include "matrix.h"
#include "network.h"
#include "newick.h"
#include "character.h"
#include "optimize.h"
#include "dumb.h"

#include <iostream>
#include <vector>
#include <set>

namespace bruteforce
{
	namespace
	{
		struct edge
		{
			network::idx_t n0, n1;
		};
		
		network::data *_network;
		std::vector<edge> _edges;
		long long _visited = 0;
		character::distance_t _best_distance;
		std::set<std::string> _best_network;
	}
	
	void next_taxon(int which)
	{
		if (which == _network->mtx_taxons)
		{	
			unsigned int d = optimize::optimize(_network);
			if (d < _best_distance)
			{
				std::cout << " -> Best new distance " << d << std::endl;
				_best_network.clear();
				_best_distance = d;
				newick::print(_network);
			}
		
			if (d == _best_distance)
				_best_network.insert(newick::from_network(_network, 0));
			
			if ((_visited++ % 10000000 == 0) && _visited > 1)
				std::cout << "...visited " << (_visited/1000000) << "M networks" << std::endl;
				
			return;
		}
		
		for (unsigned int i=0;i<_edges.size();i++)
		{
			edge tmp = _edges[i];
			network::idx_t neu = network::insert(_network, _edges[i].n0, _edges[i].n1, which);
			
			_edges[i].n0 = tmp.n0;
			_edges[i].n1 = neu;
			
			edge e1;
			e1.n0 = which;
			e1.n1 = neu;
			
			edge e2;
			e2.n0 = neu;
			e2.n1 = tmp.n1;
			
			_edges.push_back(e1);
			_edges.push_back(e2);
			
			unsigned int dist = optimize::optimize(_network);
			if (dist <= _best_distance)
				next_taxon(which+1);
			
			_edges.pop_back();
			_edges.pop_back();

			// restore			
			_edges[i] = tmp;
			
			network::disconnect(_network, which);
		}		
	}

	network::data* run(matrix::data *matrix, int bound)
	{		
		_network = network::alloc(matrix);
		_visited = 0;
		_best_distance = 6553500;
		_best_network.clear();
		
		// make a few random for some kind of bound
		
		network::init(_network, 0, 1);
		
		for (int i=0;i<16;i++)
		{
			dumb::make_inplace(_network);
			unsigned int dist = optimize::optimize(_network);
			if (i == 0 || dist < _best_distance)
			{
				_best_network.clear();
				_best_network.insert(newick::from_network(_network, 0));
				_best_distance = optimize::optimize(_network);
			}
		}
		
		if (bound != -1)
			_best_distance = bound;
		
		std::cout << "Starting brute force search (branch & bound) with start dist = " << _best_distance << std::endl;
		network::init(_network, 0, 1);
		
		_edges.clear();

		edge e;
		e.n0 = 0;
		e.n1 = 1;
		_edges.push_back(e);
		
		next_taxon(2);
		
		network::free(_network);
		
		std::cout << "Enumerated " << _visited << " networks." << std::endl;
		
		if (_best_distance < bound || bound == -1)
		{
			std::cout << _best_network.size() << " networks share minimal distance " << _best_distance << std::endl;
			for (std::set<std::string>::iterator i=_best_network.begin();i!=_best_network.end();i++)
			{
				std::cout << (*i) << std::endl;
			}
		}
		else
		{
			std::cout << "Found no networks with d < " << bound << std::endl;
		}
		
		return 0;
	}
}

