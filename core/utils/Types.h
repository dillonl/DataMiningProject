#ifndef CORE_UTILS_TYPES_H
#define CORE_UTILS_TYPES_H

#include <cstdint>
#include <limits>

namespace dmp
{
	//typedef uint64_t InternalKmer;
	typedef uint8_t InternalKmer;
	typedef uint32_t position;

	struct uint256_t
	{
		uint64_t bits[4];
	};

//#define KMER_SIZE 25ULL
#define KMER_SIZE 4ULL
#define KMER_SHIFTER_SIZE (KMER_SIZE * 2)
#define KMER_COUNT_INC 1ULL
#define KMER_MASK std::numeric_limits< InternalKmer >::max() >> ((sizeof(InternalKmer) * 8) - KMER_SHIFTER_SIZE)
}

#endif //CORE_UTILS_TYPES_H
