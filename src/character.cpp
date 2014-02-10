#include "character.h"
#include <string>
#include <iostream>

namespace character
{
	struct buf
	{
		state_t *ptr;
		unsigned int pitch;
	};

	state_t interpret(char t)
	{
		if (t >= '0' && t <= '9')
			return t - '0';
		return UNKNOWN_CHAR_VALUE;
	}
	
	char uninterpret(state_t t)
	{
		if (t == UNKNOWN_CHAR_VALUE)
			return '?';
		else return '0' + t;
	}
	
	buf* alloc(unsigned int characters, unsigned int buffers, state_t **out)
	{
		buf* b = new buf();
		b->pitch = characters;
		b->ptr = new state_t[buffers * b->pitch];
		
		for (unsigned int i=0;i<buffers;i++)
			out[i] = b->ptr + i * b->pitch;
			
		return b;
	}
	
	void free(buf *b)
	{
		delete [] b->ptr;
		delete b;
	}

	void from_string(const char *str, state_t *out)
	{
		for (int i=0;i<strlen(str);i++)
			out[i] = interpret(str[i]);
	}
	
	const char *to_string(state_t *in, int characters)
	{
		static char tmp[1024];
		for (int i=0;i<characters;i++)
			tmp[i] = uninterpret(in[i]);
		tmp[characters] = 0;
		return tmp;
	}
	
	//
	distance_t distance(state_t *a, state_t *b)
	{
		return 0;
	}
	
	//
	void threesome(state_t *a, state_t *b, state_t *c, state_t *out)
	{
	
	}
}