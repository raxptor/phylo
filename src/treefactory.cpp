#include <iostream>
#include <cstdlib>
#include <cstring>
#include <set>
#include <vector>

#include "optimize.h"
#include "matrix.h"
#include "dumb.h"
#include "smart.h"
#include "network.h"
#include "newick.h"

#include <mtw/mersenne-twister.h>
#include <algorithm>

void help()
{
	std::cout << "No help for this program, read the source" << std::endl;
	exit(0);
}

void rounds(network::data *net, int count, int method)
{
	const char *mn[] = 
	{ "dumb", "smart" };
	
	int sum = 0;
	std::cout << "Running " << count << " rounds of " << mn[method] << "..." << std::endl;
	
	std::vector<int> results;

	for (int i=0;i<count;i++)
	{
		switch (method)
		{
			case 0:
				dumb::make_inplace(net);
				break;
			case 1:
				smart::make_inplace(net);
				break;
			default:
				std::cout << "no such method" << std::endl;
				return;
		}
		
		const int score = optimize::optimize(net);
		results.push_back(score);
		sum += score;
	}
	
	float avg = float(sum)/float(count);

	std::cout << " => Average tree length=" << avg << " median=" << results[results.size()/2] << std::endl;
}

int main(int argc, const char **argv)
{
	const char *method = 0;
	const char *print = 0;
	
	mtw_seed(0);
	
	if (argc < 2)
		help();
	
	character::init();
	optimize::init();

	matrix::data *mtx = matrix::load(argv[argc-1]);
	if (!mtx) 
	{
		std::cerr << "Matrix failed to load. Aborting." << std::endl;
		return -1;
	}
	
	network::data *net = network::alloc(mtx);

	const int count = 2000;
	rounds(net, count, 0);
	rounds(net, count, 1);

	return 0;
}
