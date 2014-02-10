
#include "matrix.h"
#include "tree.h"
#include "newick.h"

#include <iostream>
#include <vector>
#include <set>

namespace bruteforce
{
	namespace
	{
		tree::data *_tree;
		std::vector<tree::idx_t> _edges;
		std::vector<tree::idx_t> _left;
		long long _visited = 0;
		bool _dontpick[1000];
	}
	
	void next_taxon();
	
	void insert(tree::idx_t taxon)
	{
		for (unsigned int i=0;i<_edges.size();i++)
		{
			if (_edges[i] < taxon)
				continue;
				
			tree::idx_t neu = tree::insert(_tree, _edges[i], taxon);

			_edges.push_back(neu);
			_edges.push_back(taxon);
			next_taxon();
			
			_edges.pop_back();
			_edges.pop_back();
			tree::disconnect(_tree, taxon);
		}
	}
	
	void next_taxon()
	{
		if (_left.empty())
		{
			if ((_visited++ % 10000000 == 0) && _visited > 1)
				std::cout << "...visited " << (_visited/1000000) << "M trees" << std::endl;
				
			return;
		}
		
		//
		for (int i=0;i<_left.size();i++)
		{
			const tree::idx_t which = _left[i];

			_left[i] = _left.back();
			_left.pop_back();
			
			insert(which);
			
			_left.push_back(_left[i]);
			_left[i] = which;
		}
	}

	tree::data* run(matrix::data *matrix)
	{		
		_tree = tree::alloc(matrix);
		_visited = 0;
		
		// for now always keep the 0 as root
		for (int i=1;i<matrix->taxons;i++)
		{
			_left.clear();
			_edges.clear();
			_edges.push_back(i);
			
			tree::init(_tree, 0, i);
			std::cout << "Matrix with next-pick [" << i << "] visited=" << _visited << std::endl;
			
			for (int j=1;j<matrix->taxons;j++)
			{
				if (j != i)
					_left.push_back(j);
			}
			next_taxon();
		}
		
		tree::free(_tree);
		
		std::cout << "Enumerated " << _visited << " trees." << std::endl;

/*		
		std::set<std::string>::iterator i = _trees.begin();
		while (i != _trees.end())
			std::cout << *i++ << std::endl;
*/
		return 0;
	}
}

