#include "character.h"
#include <string>
#include <iostream>
#include <cstring>

namespace character
{
	struct buf
	{
		state_t *ptr;
		unsigned int pitch;
		unsigned int size;
	};
	
	static const state_t res[8] = {
		0, 1, 2, 1,
		0, 0, 0, 0
	};
	
	static distance_t res8[256];
	
	void init()
	{
		for (int i=0;i<8;i++)
			for (int j=0;j<8;j++)
				res8[i*16+j] = res[i] + res[j];
	}

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
		b->pitch = 4*((characters + 5)/4);
		b->ptr = new state_t[buffers * b->pitch];
		b->size = buffers * b->pitch * sizeof(state_t);
		
		memset(b->ptr, 0x00, sizeof(state_t) * buffers * b->pitch);
		
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
		for (int i=0;i<characters;i+=2)
		{
			const int _a = a[i] + (a[i+1] << 4);
			const int _b = b[i] + (b[i+1] << 4);
			sum += res8[_a^_b];
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
	
	//
	void toggle_boost(state_t *target, int *picks, int count)
	{
		// TODO: Fix this, add internal multiplier or whatever
		for (int i=0;i<count;i++)
		{
			int w = picks[i];
			if (target[w] == 1)
				target[w] = 2;
			else if (target[w] == 2)
				target[w] = 1;
		}
	}
	
}