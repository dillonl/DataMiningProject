#ifndef CORE_CONTAINERS_MINTOMICKMERSET_HPP
#define CORE_CONTAINERS_MINTOMICKMERSET_HPP

#include "KmerHash.hpp"
#include "IKmerSet.h"

#include <memory>
#include <unordered_set>

#include <mintomic/mintomic.h>

namespace dmp
{
	class MintomicKmerSet : public IKmerSet
	{
	public:
        typedef std::shared_ptr< MintomicKmerSet > SharedPtr;
        MintomicKmerSet()
		{
		}

		~MintomicKmerSet()
		{
			// static std::mutex lock;
			// std::lock_guard< std::mutex > guard(lock);
			// std::cout << "set destroyed" << std::endl;
			delete[] m_kmer_set;
			m_total_set_size = 0;
			mint_store_64_relaxed(&m_current_set_size, 0);
		}

		/*
		//----------------------------------------------
        // from code.google.com/p/smhasher/wiki/MurmurHash3
		inline static uint64_t integerHash(uint64_t h)
		{
			h ^= h >> 33;
			h *= 0xff51afd7ed558ccd;
			h ^= h >> 33;
			h *= 0xc4ceb9fe1a85ec53;
			h ^= h >> 33;
			return h;
		}
		*/
#define BIG_CONSTANT(x) (x##LLU)
		static inline uint64_t integerHash( uint64_t k )
		{
			k ^= k >> 33;
			k *= BIG_CONSTANT(0xff51afd7ed558ccd);
			k ^= k >> 33;
			k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
			k ^= k >> 33;

			return k;
		}

		void addKmer(InternalKmer key) override
		{
			InternalKmer value = 1;
			if (key == 0) { return; }
			assert(key != 0);
			assert(value != 0);

			for (uint32_t idx = integerHash(key);; idx++)
			{
				idx &= m_total_set_size - 1;

				uint32_t prevKey = mint_compare_exchange_strong_32_relaxed(&m_entries[idx].key, 0, key);
				if ((prevKey == 0) || (prevKey == key))
				{
					mint_store_32_relaxed(&m_entries[idx].value, value);
					return;
				}
			}
		}

		uint64_t getKmerCount(InternalKmer internalKmer) override
		{
			/*
			auto iter = m_kmer_set.find(internalKmer);
			return (iter != m_kmer_set.end()) ? ((*iter) >> KMER_SHIFTER_SIZE) : 0;
			*/
			return 0;
		}

		size_t getSetSize() override
		{
			return mint_load_64_relaxed(&m_current_set_size);
		}

		void addAllKmersToPassedInSet(IKmerSet::SharedPtr kmerSetPtr) override
		{
		}

		void resize(uint64_t arraySize) override
		{
			assert((arraySize & (arraySize - 1)) == 0);   // Must be a power of 2
			m_arraySize = arraySize;
			m_entries = new Entry[arraySize];
			memset(m_entries, 0, sizeof(Entry) * m_arraySize);

			/*
			m_kmer_set = new mint_atomic64_t[size];
		    memset(m_kmer_set, 0, sizeof(mint_atomic64_t) * size);
			m_total_set_size = size;
			mint_store_64_relaxed(&m_current_set_size, 0);
			*/
			// m_current_set_size = size;
		}

	private:
		struct Entry
		{
			mint_atomic32_t key;
			mint_atomic32_t value;
		};

		Entry* m_entries;
		uint32_t m_arraySize;

		// std::unordered_set< InternalKmer, KmerHash, KmerKeyEqual > m_kmer_set;
		mint_atomic64_t* m_kmer_set;
		size_t m_total_set_size;
		mint_atomic64_t m_current_set_size;

	};
}

#endif // CORE_CONTAINERS_MINTOMICKMERSET_HPP
