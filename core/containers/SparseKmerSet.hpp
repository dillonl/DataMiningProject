#ifndef CORE_CONTAINERS_SPARSEKMERSET_HPP
#define CORE_CONTAINERS_SPARSEKMERSET_HPP

#include "IKmerSet.h"
#include "KmerHash.hpp"

#include "utils/Types.h"

#include <memory>

#include <google/sparse_hash_set>

#include <boost/noncopyable.hpp>

namespace dmp
{

	class SparseKmerSet : public IKmerSet
	{
	public:
		typedef std::shared_ptr< SparseKmerSet > SharedPtr;
		SparseKmerSet() : m_bit_mask(KMER_MASK)
		{
		}

		~SparseKmerSet()
		{
			m_kmer_set.clear();
		}

		void addKmer(InternalKmer internalKmer) override
		{
			auto iter = m_kmer_set.find(internalKmer);
			if (iter == m_kmer_set.end())
			{
				internalKmer = (internalKmer & m_bit_mask);
				m_kmer_set.insert(internalKmer);
				iter = m_kmer_set.find(internalKmer);
			}
			auto iterPtr = const_cast< InternalKmer* >(&(*iter));
			*iterPtr = (internalKmer | ((((*iter) >> KMER_SHIFTER_SIZE) + KMER_COUNT_INC) << KMER_SHIFTER_SIZE));
		}

		uint64_t getKmerCount(InternalKmer internalKmer) override
		{
			auto iter = m_kmer_set.find(internalKmer);
			return (iter != m_kmer_set.end()) ? ((*iter) >> KMER_SHIFTER_SIZE) : 0;
		}

		size_t getSetSize() override
		{
			return m_kmer_set.size();
		}

		void addAllKmersToPassedInSet(IKmerSet::SharedPtr kmerSetPtr) override
		{
		}

		void resize(uint64_t size) override
		{
			m_kmer_set.resize(size);
		}

	private:
		google::sparse_hash_set< InternalKmer, KmerHash, KmerKeyEqual > m_kmer_set;
		InternalKmer m_bit_mask;

	};
}

#endif // CORE_CONTAINERS_SPARSEKMERSET_HPP
