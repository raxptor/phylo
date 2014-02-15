#include "matrix.h"
#include "character.h"

#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>
#include <cstring>

namespace matrix
{
	struct idata
	{
		std::vector<std::string> c_name, t_name;
		std::vector<std::string> taxon_lines;
		int maxtlen;
	};
	
	data *alloc()
	{
		data *d = new data();

		d->characters = 0;
		d->taxons = 0;
		
		memset(&d->ordered, 0x00, sizeof(cgroup));
		memset(&d->unordered, 0x00, sizeof(cgroup));
		
		d->_id = new idata();
		d->_id->maxtlen = 12;
		return d;
	}
	
	void free(data *d)
	{
		delete [] d->ordered.submatrix;
		delete [] d->unordered.submatrix;
		delete [] d->ordered.weights;
		delete [] d->unordered.weights;
		delete d->_id;
		delete d;
	}
	
	void sort_characters(character::state_t *mtx, data *out)
	{
		// rows in mtx are taxons
		
		std::vector<int> ordered;
		std::vector<int> unordered;

		out->ordered.bits = 0;
		out->unordered.bits = 0;
		
		int discarded = 0;
		
		for (int i=0;i<out->characters;i++)
		{
			bool used[256];
			int count = 0;
			memset(used, 0x00, 256 * sizeof(bool));

			// make the range 0-based		
			int min = 100000;
			for (int j=0;j<out->taxons;j++)
			{
				character::state_t val = mtx[j * out->characters + i];
				if (val == character::UNKNOWN_CHAR_VALUE)
					continue;
					
				if (val < min) 
					min = val;
			}
			
			int bits = 0;
			for (int j=0;j<out->taxons;j++)
			{
				character::state_t & val = mtx[j * out->characters + i];
				if (val == character::UNKNOWN_CHAR_VALUE)
					continue;
				
				val -= min;
				
				if ((val+1) > bits)
					bits = (val+1);
				
				if (!used[val])
				{
					used[val] = true;
					count++;
				}
			}
				
			if (count == 1)
			{
				discarded++;
				std::cout << "[matrix] - discarding character " << i << " because it has only 1 value" << std::endl;
				continue;
			}
			
			cgroup *tgt;
			
//			if (used[0] && used[1] && (count == 2 || (count == 3 && used[character::UNKNOWN_CHAR_VALUE])))
//			{
				tgt = &out->unordered;
				unordered.push_back(i);
//			}
//			else
//			{
//				tgt = &out->ordered;
//				ordered.push_back(i);
//			}
			
			if (bits > tgt->bits)
				tgt->bits = bits;
		}
		
		out->ordered.count = ordered.size();
		out->unordered.count = unordered.size();

		out->ordered.submatrix = new character::state_t[ordered.size() * out->taxons];
		out->unordered.submatrix = new character::state_t[unordered.size() * out->taxons];
		out->ordered.weights = new int[ordered.size()];
		out->unordered.weights = new int[unordered.size()];
		
		// Rotate all characters into their matrix
		//
		// c   taxons
		// 0: [abcdefg] 
		// 1: [abcdefg]
		//
		for (unsigned int i=0;i<ordered.size();i++)
		{
			for (unsigned int j=0;j<out->taxons;j++)
				out->ordered.submatrix[i * out->taxons + j] = mtx[j * out->characters + ordered[i]];
				
			out->ordered.weights[i] = 1;
		}

		for (unsigned int i=0;i<unordered.size();i++)
		{
			for (unsigned int j=0;j<out->taxons;j++)
				out->unordered.submatrix[i * out->taxons + j] = mtx[j * out->characters + unordered[i]];
				
			out->unordered.weights[i] = 1;
		}
		
		std::cout << "[matrix] - Processed " << out->characters << " characters into " << out->ordered.count << " ordered and " << out->unordered.count << " unordered, " << discarded  << " discarded." << std::endl;
	}
	
	data* load(const char *fn)
	{
		std::ifstream f(fn);
		if (!f.good())
			return 0;
			
		std::vector<std::string> t;
		std::vector<std::string> c;
		
		while (!f.eof())
		{
			std::string line;
			std::getline(f, line);
			
			int w1 = line.find_first_of(' ');
			int w2 = line.find_first_of(9);
			
			if (w1 != std::string::npos && w2 != std::string::npos && w2 < w1)
				w1 = w2;
			if (w1 == std::string::npos)
				w1 = w2;
			
			// no space or tab on this line
			if (w1 == std::string::npos)
				continue;
			
			t.push_back(line.substr(0, w1));
			
			std::string tmp;
			
			// filter out characters
			for (int i=w1+1;i<line.size();i++)
			{
				if ((line[i] >= '0' && line[i] <= '9') || (line[i] == '?' || line[i] == '-'))
					tmp.push_back(line[i]);
			}
			
			c.push_back(tmp);
		}
		
		if (t.size() < 2)
		{
			std::cerr << "Error: Matrix contains less than 2 taxons" << std::endl;
			return 0;
		}
		
		for (int i=1;i<t.size();i++)
		{
			if (c[i].size() != c[0].size())
			{
				std::cerr << "Error: " << t[i] << " has " << c[i].size() << " characters, but " << t[0] << " has " << c[0].size() << ". Make it a matrix and try again." << std::endl;
				return 0;
			}
		}
		
		
		data *d = alloc();
		d->_id->maxtlen = 0;
		d->characters = c[0].size();
		d->taxons = t.size();
		
		// Temporary buffer for preprocessing
		character::state_t *mtx = new character::state_t[d->characters * d->taxons];
		
		for (unsigned int i=0;i<t.size();i++)
		{
			// max taxon name length for pretty printing
			if (t[i].size() > d->_id->maxtlen)
				d->_id->maxtlen = t[i].size();
			
			// --
			character::from_string(c[i].c_str(), &mtx[i * d->characters]);
		}
		
		sort_characters(mtx, d);
		
		delete [] mtx;
		
		d->_id->t_name = t;

		std::cout << "[matrix] - Done loading matrix with " << t.size() << " taxons and " << c[0].size() << " characters." << std::endl;
		return d;
	}
	
	const char *character_name(data *d, unsigned int index)
	{
		if (index < d->_id->c_name.size())
			return d->_id->c_name[index].c_str();
		else
			return "no-name";
	}

	const char *taxon_name(data *d, unsigned int index)
	{
		if (index < d->_id->t_name.size())
			return d->_id->t_name[index].c_str();
		else
			return "no-name";
	}

	void print(data *d)
	{
	
	}
}
