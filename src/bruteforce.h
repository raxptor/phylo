#ifndef __BRUTEFORCE_H__
#define __BRUTEFORCE_H__

namespace tree { struct data; }
namespace matrix { struct data; }

namespace bruteforce
{
	// brute forces the maximally parsimonious tree
	tree::data* run(matrix::data *matrix);
}

#endif
