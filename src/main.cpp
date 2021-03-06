#include <iostream>
#include <cstdlib>
#include <cstring>
#include <set>

#include "optimize.h"
#include "matrix.h"
#include "dumb.h"
#include "smart.h"
#include "bruteforce.h"
#include "newick.h"
#include "tbr.h"
#include "ratchet.h"
#include "network.h"

#include <mtw/mersenne-twister.h>
#include <fstream>

extern long long g_treeTiming;

void help()
{
	std::cout << "Run: phylo [args] matrix-file" << std::endl; 
	std::cout << "--seed <seed>  - Specify random seed" << std::endl;
	std::cout << "--run <method> - Run analysis with method [bruteforce, ratchet]" << std::endl;
	std::cout << "--print <what> - Print [matrix]" << std::endl;
	std::cout << "--timing <MTrees> - stop after specified number million trees" << std::endl;
	std::cout << std::endl;
	
	exit(0);
}

namespace tbr { extern long long networks; }

void print_result(tbr::output *out)
{
	std::cout << std::endl;
	std::cout << "Done. I have " << out->equal_length.size() << " networks of length " << out->length << std::endl;
	for (int i=0;i<out->newick.size()&&i<200;i++)
		std::cout << out->newick[i] << std::endl;
		
	std::cout << "Networks searched: " << (tbr::networks/1000000) << "M" << std::endl;
	std::cout << "Printed networks are of length " << out->length << std::endl;
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
	int rounds = 10;
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
		else if (!strcmp(argv[c], "--rounds"))
		{
			rounds = atoi(argv[c+1]);
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
		else if (!strcmp(argv[c], "--timing"))
		{
			g_treeTiming = atoi(argv[c+1])*1000000;
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
	optimize::init();
	
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
		
	if (method && !strcmp(method, "optimize"))
	{
		network::data *tmp = dumb::make(mtx);
		optimize::optimize(tmp);
		network::print_characters(tmp);
	}

	if (method && !strcmp(method, "smart"))
	{
		network::data *tmp = smart::make(mtx);
		std::cout << "smart length = " << optimize::optimize(tmp) << std::endl;
		network::print_characters(tmp);
	}
	
	if (method && !strcmp(method, "tbr"))
	{
		network::data *tmp = smart::make(mtx);
		
		tbr::output output;
		output.best_network = network::alloc(mtx);
		network::copy(output.best_network, tmp);
		output.length = optimize::optimize(output.best_network);
		
		for (int i=0;i<rounds;i++)
		{
			smart::make_inplace(tmp);
			while (tbr::run(tmp, &output))
			{
				network::copy(tmp, output.best_network);
			}
		}
		
		network::copy(output.best_network, tmp);
		
		print_result(&output);
	}

	if (method && !strcmp(method, "ratchet"))
	{
		const int num = 5; // 1 best and 1 speculative (best will try ratchet differently, speculative will wander around)
		network::data* nw[num];
		for (int i=0;i<num;i++)
			nw[i] = smart::make(mtx);

		tbr::output out;
		out.best_network = network::alloc(mtx);
		network::copy(out.best_network, nw[0]);
		out.length = optimize::optimize(out.best_network);
		
		std::cout << "Starting ratchet with random tree dist=" << out.length << std::endl;

		int maxnew = 10;
		
		for (int i=0;i<rounds/num;i++)
		{
			int oc = out.equal_length.size();
			
			unsigned int old = out.length;
			
			for (int j=0;j<num;j++)
			{
				for (int k=0;k<10;k++)
				{
					if (tbr::run(nw[j], &out))
					{
						// keep the best in j=0 and let the others be speculative
						if (j == 0)
							network::copy(nw[j], out.best_network);
					}
					else
					{
						break;
					}
				}
			}

			if (out.length < old)
			{
				std::cout << "[tbr pre-ratchet] - i have " << out.equal_length.size() << " of the same length (" << out.length << ")" << std::endl;
				std::cout << "[tbr pre-ratchet] - here is one: " << *(out.equal_length.begin()) << std::endl;
			}

			for (int j=0;j<num;j++)
			{
				ratchet::run(nw[j], &out);
			}

			// eliminate dupes						
			char tmpbuf[num][16000];
			for (int i=0;i<num;i++)
			{
				char *p = &tmpbuf[i][0];
				network::sort(nw[i]);
				network::to_string_2(nw[i], &p, 0, network::NOT_IN_NETWORK);
				*p = 0;
			}
				
			for (int i=0;i<num;i++)
			{
				for (int j=i+1;j<num;j++)
				{
					if (!strcmp(tmpbuf[i], tmpbuf[j]))
					{
						if (maxnew-- > 0)
							smart::make_inplace(nw[j]);
					}
				}
			}
					

			if (oc != out.equal_length.size())
			{
				std::cout << "[ratchet] - i have " << out.equal_length.size() << " of the same length (" << out.length << ")" << std::endl;
			}
		}
		print_result(&out);
	}
	
	matrix::free(mtx);	
	return 0;
}
