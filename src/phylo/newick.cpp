#include <iostream>
#include <cstdlib>
#include <algorithm>

#include "network.h"

namespace newick
{
	struct subparts 
	{
		const char *data;
		int min_tax;
		inline bool operator<(subparts const &b) const
		{
			return min_tax < b.min_tax;
		}
	};
	
	std::string subnet(network::data *d, int where, int from, int *min_tax, bool isroot = true)
	{
		network::idx_t c[3];
		c[0] = d->network[where].c0;
		c[1] = d->network[where].c1;
		c[2] = d->network[where].c2;
				
		subparts sp[3];
		std::string st[3];
		int out = 0;
		
		*min_tax = where;
		
		for (int i=0;i<3;i++)
		{			
			if (c[i] != from && c[i] != network::NOT_IN_NETWORK)
			{
				st[out] = subnet(d, c[i], where, &sp[out].min_tax, false);
				sp[out].data = st[out].c_str();
				if (sp[out].min_tax < *min_tax)
					*min_tax = sp[out].min_tax;
				++out;
			}
		}
		
		std::sort(&sp[0], &sp[out]);
		
		std::string cont;
		bool first = true;
		
		for (int i=0;i<out;i++)
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
			cont.append(sp[i].data);
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
		int tmp;
		return std::string("(") + subnet(d, root, network::NOT_IN_NETWORK, &tmp, true) + ");";
	}

	void print(network::data *d, int root)
	{
		std::cout << from_network(d, root) << std::endl;
	}
}
