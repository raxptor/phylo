#ifndef __TBR_H__
#define __TBR_H__

namespace network { struct data; }

namespace tbr
{
	struct output
	{
		network::data *best_network;
	};

	int run(network::data *net, output *out);
}

#endif
