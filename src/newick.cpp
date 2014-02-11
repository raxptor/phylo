#include <iostream>
#include <cstdlib>

#include "network.h"
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

		return cont;		
	}
	
	std::string subnet(network::data *d, int where, int from, bool isroot = true)
	{
		network::idx_t c[3];
		c[0] = d->network[where].c0;
		c[1] = d->network[where].c1;
		c[2] = d->network[where].c2;

		std::string cont;
		bool first = true;
		char tmp[64];

		for (int i=0;i<3;i++)
		{
			if (c[i] != from && c[i] != network::NOT_IN_NETWORK)
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
				
				cont.append(subnet(d, c[i], where, false));
			}
		}

		if (!first)
			cont.append(")");
		
		if (isroot)
			cont.append(",");
		
		if (where < d->matrix->taxons)
		{
			cont.append(matrix::taxon_name(d->matrix, where));
		}
		
		return cont;
	}
	
	
	std::string from_network(network::data *d, int root)
	{
		return std::string("(") + subnet(d, root, network::NOT_IN_NETWORK, true) + ");";
	}

	std::string from_tree(tree::data *d)
	{
		return std::string("(") + subtree(d, 0) + ");";
	}

	void print(network::data *d, int root)
	{
		std::cout << from_network(d, root) << std::endl;
	}

	void print(tree::data *d)
	{
		std::cout << "Newick (tree): " << from_tree(d) << std::endl;
	}
}
