#ifndef __RATCHET_H__
#define __RATCHET_H__

namespace network { struct data; }
namespace tbr { struct output; }

namespace ratchet
{
	void run(network::data *data, tbr::output * out);
}

#endif
