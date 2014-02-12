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

#include <mtw/mersenne-twister.h>

int main(int argc, const char **argv)
{
	
	if (argc < 2)
	{
		std::cerr << "Error: no matrix specified." << std::endl;
		return -1;
	}
	
	if (argc > 2)
		seed(atoi(argv[2]));
	
	character::init();
	
	matrix::data *mtx = matrix::load(argv[1]);
	if (!mtx) 
	{
		std::cerr << "Matrix failed to load. Aborting." << std::endl;
		return -1;
	}
	
	matrix::print(mtx);

	const int num = 10;
	network::data* nw[num];
	for (int i=0;i<num;i++)
		nw[i] = dumb::make(mtx);

	tbr::output out;
	out.best_network = network::alloc(mtx);
	network::copy(out.best_network, nw[0]);

	int rot = 0;
	for (int i=0;i<10000;i++)
	{
		for (int j=0;j<num;j++)
			ratchet::run(nw[j], &out);
			
		rot = (rot + 1) % num;
		if (nw[rot]->dist >= out.best_network->dist)
			dumb::make_inplace(nw[rot]);
	}

	std::cout << std::endl;
	std::cout << "Done. Best network (" << out.best_network->dist << ") ==> ";
	newick::print(out.best_network);
/*
	std::cout << "Bruteforcing the solution..." << std::endl;
	bruteforce::run(mtx);
*/
	matrix::free(mtx);
	
	return 0;
}
