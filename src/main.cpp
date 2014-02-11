#include <iostream>
#include <cstdlib>

#include "matrix.h"
#include "dumb.h"
#include "tree.h"
#include "bruteforce.h"
#include "newick.h"
#include "tbr.h"
#include "ratchet.h"
#include "network.h"

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

	network::data *nw = dumb::make(mtx);

	std::cout << "Network starting dist: " << nw->dist << " => ";	
	newick::print(nw);

//	for (int i=0;i<2000;i++)
//		tbr::run(nw);

	for (int i=0;i<2000;i++)
		ratchet::run(nw);

	std::cout << std::endl;
	std::cout << "Done. Best network (" << nw->dist << ") ==> ";
	newick::print(nw);
/*
	std::cout << "Bruteforcing the solution..." << std::endl;
	bruteforce::run(mtx);
*/
	matrix::free(mtx);
	
	return 0;
}
