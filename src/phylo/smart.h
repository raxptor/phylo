#ifndef __SMART_H__
#define __SMART_H__

namespace network { struct data; }
namespace matrix { struct data; }

namespace smart
{
	// creates a completely random network
	void make_inplace(network::data *network);
	
	// allocates new network for the matrix and returns one
	network::data *make(matrix::data *mtx);
};

#endif
