#include "dumb.h"
#include "tree.h"
#include "matrix.h"

#include <vector>
#include <cstdlib>
#include <iostream>

namespace dumb
{
	tree::data *make(matrix::data *mtx)
	{
		if (mtx->characters < 2)
		{
			std::cerr << "Trying to make a tree with less than 2 characters!" << std::endl;
			return 0;
		}
		
		tree::data *t = tree::alloc(mtx);
		make_inplace(t);
		return t;
	}
	
	void make_inplace(tree::data *t)
	{
		tree::init(t, 0, 1);
		std::vector<tree::idx_t> edges;
		
		edges.reserve(t->matrix->taxons * 2);
		edges.push_back(1); 
		
		for (int i=2;i<t->matrix->taxons;i++)
		{
			const int pick = rand() % edges.size();
			edges.push_back(tree::insert(t, edges[pick], i));
			edges.push_back(i);
		}
	}
}
