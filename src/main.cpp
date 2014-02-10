#include <iostream>
#include <cstdlib>

#include "matrix.h"
#include "dumb.h"
#include "tree.h"

int main(int argc, const char **argv)
{
	if (argc < 2)
	{
		std::cerr << "Error: no matrix specified." << std::endl;
		return -1;
	}
	
	if (argc > 2)
		srand(atoi(argv[2]));
	
	matrix::data *mtx = matrix::load(argv[1]);
	if (!mtx) 
	{
		std::cerr << "Matrix failed to load. Aborting." << std::endl;
		return -1;
	}
	
	tree::data *tr = dumb::make(mtx);
	
	if (tr)
	{
		tree::print_newick(tr);
	}
	
	tree::free(tr);
	matrix::free(mtx);
	
	return 0;
}
