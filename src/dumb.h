#ifndef __DUMB_H__
#define __DUMB_H__

namespace tree { struct data; }
namespace matrix { struct data; }

namespace dumb
{
	// creates a completely random tree
	void make_inplace(tree::data *tree);
	
	// allocates new tree for the matrix and returns one
	tree::data *make(matrix::data *mtx);
};

#endif
