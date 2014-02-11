
#include "matrix.h"
#include "tree.h"
#include "newick.h"
#include "character.h"

#include <iostream>
#include <vector>
#include <set>

namespace bruteforce
{
	namespace
	{
		tree::data *_tree;
		std::vector<tree::idx_t> _edges;
		long long _visited = 0;
		character::distance_t _best_distance;
		std::vector<std::string> _best_tree;
	}
	
	void next_taxon(int which);
	
	void insert()
	{

	}
	
	void next_taxon(int which)
	{
		if (which == _tree->mtx_taxons)
		{
			if (_tree->dist > _best_distance)
			{
			
			}
			else
			{
				if (_tree->dist < _best_distance)
				{
					std::cout << " -> Best new distance " << _tree->dist << std::endl;
					_best_tree.clear();
					_best_distance = _tree->dist;
				}
				_best_tree.push_back(newick::from_tree(_tree));
			}
			
			if ((_visited++ % 10000000 == 0) && _visited > 1)
				std::cout << "...visited " << (_visited/1000000) << "M trees" << std::endl;
				
			return;
		}
		
		for (unsigned int i=0;i<_edges.size();i++)
		{
			tree::idx_t neu = tree::insert(_tree, _edges[i], which);
			_edges.push_back(neu);
			_edges.push_back(which);
			
			next_taxon(which+1);
			
			_edges.pop_back();
			_edges.pop_back();
			tree::disconnect(_tree, which);
		}		
	}

	tree::data* run(matrix::data *matrix)
	{		
		_tree = tree::alloc(matrix);
		_visited = 0;
		_best_distance = 65535;
		_best_tree.clear();
		
		tree::init(_tree, 0, 1);
		
		_edges.clear();
		_edges.push_back(1);
		next_taxon(2);
		
		tree::free(_tree);
		
		std::cout << "Enumerated " << _visited << " trees." << std::endl;
		
		std::cout << _best_tree.size() << " trees share minimal distance " << _best_distance << std::endl;
		for (unsigned int i=0;i<_best_tree.size() && i < 10;i++)
		{
			std::cout << "Tree " << i << ": " << _best_tree[i] << std::endl;
		}
		
		return 0;
	}
}

