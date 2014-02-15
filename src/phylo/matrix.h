#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "character.h"

namespace matrix
{
	struct idata;
	
	struct cinfo
	{
		int weight;
	};
	
	struct cgroup
	{
		unsigned int count;		// number of characters in group
		unsigned int bits;		// bits needed per value (1 per value)
		cinfo *info;			// additional information
		character::state_t *submatrix;  // [taxon][count] sized array
		int *weights;
	};
	
	struct data
	{
		unsigned int characters;
		unsigned int taxons;
		
		cgroup unordered;
		cgroup ordered;
		
		// internal
		idata *_id;
	};
	
	data* load(const char *file);
	data *alloc(unsigned int taxons, unsigned int characters);
	
	void free(data *d);
	void print(data *d);
	
	unsigned int taxon_count(data *d);
	
	const char *char_name(data *d, unsigned int index);
	const char *taxon_name(data *d, unsigned int index);
	
}

#endif
