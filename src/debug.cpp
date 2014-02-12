#include "network.h"
#include "newick.h"

#include <iostream>

namespace debug
{
	void mk_backtrack(network::data *d)
	{
		network::data *d2 = network::alloc(d->matrix);
		network::copy(d2, d);
		std::cout << "logentry addlog[] = {" << std::endl;
		std::cout << "{" << d2->dist << ", \"" << newick::from_network(d2, 0) << "\"}, " << std::endl;
		for (int i=d2->mtx_taxons-1;i>1;i--)
		{
			int r0, r1;
			network::disconnect(d2, i);
			std::cout << "{" << d2->dist << ", \"" << newick::from_network(d2, 0) << "\"}, " << std::endl;
		}	
		std::cout << "};" << std::endl;
		network::free(d2);
	}
}



