#ifndef CORE_CONTAINERS_KMERHASH_HPP
#define CORE_CONTAINERS_KMERHASH_HPP

#include "utils/Types.h"

#include <memory>

namespace dmp
{
	class KmerHash
	{
	public:
		KmerHash() : m_bit_mask(KMER_MASK)
		{
		}

		size_t operator()(const InternalKmer& kmer) const
		{
			// return kmer & 0x0003FFFFFFFFFFFF;
			return kmer & this->m_bit_mask;
		}
	private:
		InternalKmer m_bit_mask;
	};

	class KmerKeyEqual
	{
	public:
		KmerKeyEqual() : m_bit_mask(KMER_MASK)
		{
		}


		bool operator()(const InternalKmer& a, const InternalKmer& b) const
		{
			// return (a & 0x0003FFFFFFFFFFFF) == (b & 0x0003FFFFFFFFFFFF);
			return (a & this->m_bit_mask) == (b & this->m_bit_mask);
		}

	private:
		InternalKmer m_bit_mask;
	};
}
#endif  //CORE_CONTAINERS_KMERHASH_HPP
