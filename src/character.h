#ifndef __CHARACTER_H__
#define __CHARACTER_H__

namespace character
{
	enum 
	{
		UNKNOWN_CHAR_VALUE = 127
	};

	typedef unsigned char state_t;
	typedef unsigned short distance_t;
	
	struct buf;
	
	//
	buf* alloc(unsigned int characters, unsigned int buffers, state_t **out);

	//
	void free(buf *b);

	// parse from string like 0120301230?123
	void from_string(const char *str, state_t *out);
	const char *to_string(state_t *buf, int characters);
	
	//
	distance_t distance(state_t *a, state_t *b);
	
	//
	void threesome(state_t *a, state_t *b, state_t *c, state_t *out);
}

#endif
