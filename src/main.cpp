#include <iostream>

#include "matrix.h"

int main(int argc, const char **argv)
{
	if (argc < 2)
	{
		std::cerr << "Error: no matrix specified." << std::endl;
		return -1;
	}
	
	matrix::data *mtx = matrix::load(argv[1]);
	if (!mtx) 
	{
		std::cerr << "Matrix failed to load. Aborting." << std::endl;
		return -1;
	}
	
	matrix::print(mtx);

	std::cout << "done." << std::endl;
//	matrix::free(mtx);
	return 0;
}
