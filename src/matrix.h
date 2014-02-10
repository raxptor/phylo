#ifndef __MATRIX_H__
#define __MATRIX_H__

namespace matrix
{
	typedef unsigned char character_state_t;

	enum 
	{
		UNKNOWN_CHAR_VALUE = 127
	};

	struct idata;
	
	struct data
	{
		unsigned int characters;
		unsigned int taxons;
		character_state_t **taxonbase;
		
		// internal
		idata *_id;
	};
	
	data* load(const char *file);
	void free(data *d);
	void print(data *d);
	
	unsigned int taxon_count(data *d);
	
	const char *char_name(data *d, unsigned int index);
	const char *taxon_name(data *d, unsigned int index);
}

#endif

