#include <iostream>
#include <cstdlib>
#include <cstring>

#include "optimize.h"
#include "matrix.h"
#include "dumb.h"
#include "tree.h"
#include "bruteforce.h"
#include "newick.h"
#include "tbr.h"
#include "ratchet.h"
#include "network.h"

#include <mtw/mersenne-twister.h>

void help()
{
	std::cout << "Run: phylo [args] matrix-file" << std::endl; 
	std::cout << "--seed <seed>  - Specify random seed" << std::endl;
	std::cout << "--run <method> - Run analysis with method [bruteforce, ratchet]" << std::endl;
	std::cout << "--print <what> - Print [matrix]" << std::endl;
	std::cout << std::endl;
	exit(0);
}

int main(int argc, const char **argv)
{
	const char *method = 0;
	const char *print = 0;
	
	mtw_seed(0);
	
	if (argc == 1)
		help();
	
	// arguments with arguments
	int bound = -1;
	int c;
	for (c=1;c<argc-1;c++)
	{
		if (!strcmp(argv[c], "--run"))
		{
			method = argv[c+1];
			c++;
		}
		else if (!strcmp(argv[c], "--bound"))
		{
			bound = atoi(argv[c+1]);
			c++;
		}
		else if (!strcmp(argv[c], "--seed"))
		{
			mtw_seed(atoi(argv[c+1]));
			c++;
		}
		else if (!strcmp(argv[c], "--print"))
		{
			print = argv[c+1];
			c++;
		}
	}
	
	if (!strcmp(argv[c], "--help"))
	{
		help();
	}
	
	if (c >= argc)
	{
		std::cerr << "Error: no matrix specified. Please provide a path to the matrix" << std::endl;
		std::cerr << "Use --help for help" << std::endl;
		return -1;
	}
	
	character::init();
	
	matrix::data *mtx = matrix::load(argv[c]);
	if (!mtx) 
	{
		std::cerr << "Matrix failed to load. Aborting." << std::endl;
		return -1;
	}
	
	if (print && !strcmp(print, "matrix"))
		matrix::print(mtx);
	
	if (method && !strcmp(method, "bruteforce"))
		bruteforce::run(mtx, bound);
		
	network::data *tmp = dumb::make(mtx);
	optimize::optimize(tmp);

	if (method && !strcmp(method, "ratchet"))
	{
		const int num = 4;
		network::data* nw[num];
		for (int i=0;i<num;i++)
			nw[i] = dumb::make(mtx);

		tbr::output out;
		out.best_network = network::alloc(mtx);
		network::copy(out.best_network, nw[0]);

		for (int i=0;i<10000;i++)
		{
			for (int j=0;j<num;j++)
			{
				ratchet::run(nw[j], &out);
			}


			// eliminate dupes						
			for (int i=0;i<num;i++)
				for (int j=i+1;j<num;j++)
					if (newick::from_network(nw[i], 0) == newick::from_network(nw[j], 0))
						dumb::make_inplace(nw[j]);
				
		}

		std::cout << std::endl;
		std::cout << "Done. Best network (" << out.best_network->dist << ") ==> ";
		newick::print(out.best_network);
	}
	
	matrix::free(mtx);	
	return 0;
}
