#include <iostream>
#include <cstdlib>

#include "matrix.h"
#include "dumb.h"
#include "tree.h"
#include "bruteforce.h"
#include "newick.h"

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
	
	matrix::print(mtx);

	std::cout << "Bruteforcing the solution..." << std::endl;

	tree::data *b = bruteforce::run(mtx);
	if (b)
	{
		std::cout << " best tree => ";
		newick::print(b);
		tree::free(b);
	}
	else
	{
		std::cout << "No result." << std::endl;
	}
	
	matrix::free(mtx);
	
	return 0;
}
