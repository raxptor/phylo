#include "matrix.h"
#include "character.h"

#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>

namespace matrix
{
	struct idata
	{
		character::buf *cbuf;	
		unsigned int pitch;
		std::vector<std::string> c_name, t_name;
		int maxtlen;
	};
	
	data *alloc(unsigned int taxons, unsigned int characters)
	{
		data *d = new data();
		d->_id = new idata();

		d->taxonbase = new character::state_t*[taxons];
		d->_id->cbuf = character::alloc(characters, taxons, d->taxonbase);
		d->_id->maxtlen = 12;
		
		d->taxons = taxons;
		d->characters = characters;
		
		
		return d;
	}
	
	void free(data *d)
	{
		character::free(d->_id->cbuf);
		delete d->_id;
		delete d;
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
				if ((line[i] >= '0' && line[i] <= '9') || line[i] == '?')
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
		
		std::cout << "Loaded matrix with " << t.size() << " taxons and " << c[0].size() << " characters." << std::endl;
		
		data *d = alloc(t.size(), c[0].size());
		d->_id->maxtlen = 0;
		
		for (unsigned int i=0;i<t.size();i++)
		{
			// max taxon name length for pretty printing
			if (t[i].size() > d->_id->maxtlen)
				d->_id->maxtlen = t[i].size();
				
			character::from_string(c[i].c_str(), d->taxonbase[i]);
		}
		
		d->_id->t_name = t;
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
		for (int i=0;i<d->taxons;i++)
		{
			std::string n = taxon_name(d, i);
			std::cout << n;
			for (int j=0;j<d->_id->maxtlen-n.size();j++)
				std::cout << " ";
				
			std::cout << " => " << character::to_string(d->taxonbase[i], d->characters);
			std::cout << std::endl;
		}
		
		std::stringstream tmp;
		tmp.setf(std::ios::left, std::ios::adjustfield);

		for (int i=0;i<d->taxons;i++)
		{
			for (int j=0;j<=i;j++)
			{
				tmp.width(4);
				tmp << character::distance(d->taxonbase[i], d->taxonbase[j], d->characters);
			}
			tmp << "\n";
		}
		
		std::cout << std::endl;
		std::cout << tmp.str() << std::endl;
	}
}
