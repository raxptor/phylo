#include "matrix.h"

#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>

namespace matrix
{
	struct idata
	{
		character_state_t *allocbuf;
		unsigned int pitch;
		std::vector<std::string> c_name, t_name;
		int maxtlen;
	};
	
	character_state_t interpret(char t)
	{
		if (t >= '0' && t <= '9')
			return t - '0';
		return UNKNOWN_CHAR_VALUE;
	}
	
	char uninterpret(character_state_t t)
	{
		if (t == UNKNOWN_CHAR_VALUE)
			return '?';
		else return '0' + t;
	}
	
	data *setup(unsigned int taxons, unsigned int characters)
	{
		data *d = new data();
		d->_id = new idata();
		
		d->taxons = taxons;
		d->characters = characters;
		
		d->_id->pitch = d->characters;
		d->_id->allocbuf = new character_state_t[d->taxons * d->_id->pitch];
		
		d->taxonbase = new character_state_t*[d->taxons];
		for (unsigned int i=0;i<d->taxons;i++)
			d->taxonbase[i] = d->_id->allocbuf + i * d->_id->pitch;
			
		return d;
	}
	
	void free(data *d)
	{
		delete [] d->_id->allocbuf;
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
				if (line[i] == '0' || line[i] == '1' || line[i] == '?')
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
		
		data *d = setup(t.size(), c[0].size());
		d->_id->maxtlen = 0;
		
		for (unsigned int i=0;i<t.size();i++)
		{
			// max taxon name length for pretty printing
			if (t[i].size() > d->_id->maxtlen)
				d->_id->maxtlen = t[i].size();
				
			for (unsigned int j=0;j<c[i].size();j++)
				d->taxonbase[i][j] = interpret(c[i][j]);
		}
		
		d->_id->t_name = t;
		return d;
	}
	
	const char *character_name(data *d, unsigned int index)
	{
		if (index < d->characters)
			return d->_id->c_name[index].c_str();
		else
			return "no-name";
	}

	const char *taxon_name(data *d, unsigned int index)
	{
		assert(index < d->taxons);
		return d->_id->t_name[index].c_str();
	}

	// if allocated as original data
	bool contains_taxon(data *d, character_state_t *t)
	{
		return t >= d->_id->allocbuf && t < d->_id->allocbuf + (d->characters * d->_id->pitch); 
	}

	void print(data *d)
	{
		for (int i=0;i<d->taxons;i++)
		{
			std::string n = taxon_name(d, i);
			std::cout << n;
			for (int j=0;j<d->_id->maxtlen-n.size();j++)
				std::cout << " ";
				
			std::cout << " => ";
			for (int k=0;k<d->characters;k++)
				std::cout << uninterpret(d->taxonbase[i][k]);
			std::cout << std::endl;
		}
	}
}
