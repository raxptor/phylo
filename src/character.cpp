#include "character.h"
#include <string>
#include <iostream>

namespace character
{
	struct buf
	{
		state_t *ptr;
		unsigned int pitch;
		unsigned int size;
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
		b->size = buffers * b->pitch * sizeof(state_t);
		
		for (unsigned int i=0;i<buffers;i++)
			out[i] = b->ptr + i * b->pitch;
			
		return b;
	}

	void copy(buf *dst, buf *src)
	{
		memcpy(dst->ptr, src->ptr, dst->size);
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
	
	void copy(state_t *target, state_t *source, int characters)
	{
		for (int i=0;i<characters;i++)
			target[i] = source[i];
	}
	
	//
	distance_t distance(state_t *a, state_t *b, int characters)
	{
		int sum = 0;
		for (int i=0;i<characters;i++)
		{
			if (a[i] != UNKNOWN_CHAR_VALUE && b[i] != UNKNOWN_CHAR_VALUE)
			{
				const int diff = a[i] - b[i];
				sum += diff > 0 ? diff : -diff;
			}
		}
		return sum;
	}
	
	//
	void threesome(state_t *a, state_t *b, state_t *c, state_t *out, int characters)
	{
		for (int i=0;i<characters;i++)
		{
			if (a[i] == b[i] || a[i] == c[i])
				out[i] = a[i];
			else if (b[i] == c[i])
				out[i] = b[i];
			else
				out[i] = c[i];
		}
	}
}