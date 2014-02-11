#include <iostream>
#include <cstdlib>

#include "tree.h"

namespace newick
{
	// this is O(n^2)!
	std::string subtree(tree::data *d, int where)
	{
		bool first = true;

		std::string cont;
		for (int i=0;i<d->allocnodes;i++)
		{
			if (d->tree[i].parent == where)
			{
				if (!first) 
				{
					cont.append(",");
				}
				else
				{
					cont.append("(");
					first = false;
				}
					
				cont.append(subtree(d, i));
			}
		}
		
		if (!first)
			cont.append(")");
		
		if (!where)
			cont.append(",");
		
		if (where < d->matrix->taxons)
		{
			cont.append(matrix::taxon_name(d->matrix, where));
		}
		else
		{
			char tmp[64];
			sprintf(tmp, "htu%d", where - d->matrix->taxons);
			cont.append(tmp);
		}

		return cont;		
	}

	std::string from_tree(tree::data *d)
	{
		return std::string("(") + subtree(d, 0) + ");";
	}

	void print(tree::data *d)
	{
		std::cout << "Newick: " << from_tree(d) << std::endl;
	}
}
