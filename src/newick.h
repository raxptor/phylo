#ifndef __NEWICK_H__
#define __NEWICK_H__

#include <string>

namespace tree { struct data; };

namespace newick
{
	void print(tree::data *d);
	std::string from_tree(tree::data *d);
}

#endif
