#ifndef __TBR_H__
#define __TBR_H__

#include <string>
#include <set>
#include <vector>

namespace network { struct data; }

namespace tbr
{
	struct output
	{
		network::data *best_network;
		std::set<std::string> equal_length;
		std::vector<std::string> newick;
		unsigned int length;
	};

	int run(network::data *net, output *out);
}

#endif
