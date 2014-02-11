#ifndef __NEWICK_H__
#define __NEWICK_H__

#include <string>

namespace tree { struct data; };
namespace network { struct data; };

namespace newick
{
	void print(tree::data *d);
	void print(network::data *d, int root=0);
	
	std::string from_tree(tree::data *d);
	std::string from_network(network::data *d, int root);
}

#endif
