#ifndef DMP_CORE_PARTIALHASH_HPP
#define DMP_CORE_PARTIALHASH_HPP

#include "utils/Types.h"

namespace dmp
{
	class DMPHash
	{
	public:
		size_t operator()(const InternalKmer& kmer) const
		{
			// return kmer & 0x03FFFFFF;
			// return kmer & 0x0004000000000000;
			std::cout << "hash: " << std::bitset< 64 >(kmer & 0x0003FFFFFFFFFFFF) << std::endl;
			return kmer & 0x0003FFFFFFFFFFFF;
		}
	};
}

#endif //DMP_CORE_PARTIALHASH_HPP
