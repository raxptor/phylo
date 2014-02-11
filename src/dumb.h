#ifndef __DUMB_H__
#define __DUMB_H__

namespace network { struct data; }
namespace matrix { struct data; }

namespace dumb
{
	// creates a completely random network
	void make_inplace(network::data *network);
	
	// allocates new network for the matrix and returns one
	network::data *make(matrix::data *mtx);
};

#endif
