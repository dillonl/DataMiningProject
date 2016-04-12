#include "KmerLookup.h"
#include "parsers/AlignmentParser.hpp"

#include <iostream>

namespace dmp
{
	KmerLookup* KmerLookup::s_instance = NULL;

	KmerLookup::KmerLookup()
	{
	}

	KmerLookup::~KmerLookup()
	{
	}

	void KmerLookup::init()
	{
		std::lock_guard< std::mutex > lg(this->m_map_init_mutex);
		std::vector< std::string > alpha = { "A", "T", "C", "G" };
		uint32_t counter = 0;
		for (uint32_t i = 0; i < 4; ++i)
		{
			for (uint32_t k = 0; k < 4; ++k)
			{
				for (uint32_t j = 0; j < 4; ++j)
				{
					for (uint32_t l = 0; l < 4; ++l)
					{
						std::string kmer = alpha[i] + alpha[k] + alpha[j] + alpha[l];
						InternalKmer iKmer = AlignmentParser::ParseKmer(kmer);
						m_order_mapping.emplace(iKmer, ++counter);
					}
				}
			}
		}
	}

	std::vector< InternalKmer > KmerLookup::getOptimalKmerSubset(const std::vector< InternalKmer >& kmers)
	{
		// static std::mutex lock;
		// std::lock_guard< std::mutex > l(lock);
		std::vector< InternalKmer > optimalKmers;
		for (auto i = 0; i < kmers.size(); i += 4)
		{
			optimalKmers.push_back(kmers[i]);
		}
		return optimalKmers;
	}

	std::vector< InternalKmer > KmerLookup::getOptimalKmerSubsetOrdered(const std::vector< InternalKmer >& kmers)
	{
		// static std::mutex l;
		// std::lock_guard< std::mutex > lock(l);
		std::vector< InternalKmer > optimalKmers;
		auto internalKmerSize = 4;

		auto orderNumber = 255;
		size_t orderIdx = 0;
		for (size_t i = 0; i < kmers.size(); ++i)
		{
			std::string tmpKmer = AlignmentParser::UnParseKmer(kmers[i]);
			auto orderNumberTmp = m_order_mapping[kmers[i]];
			if (orderNumberTmp < orderNumber)
			{
				orderNumber = orderNumberTmp;
				orderIdx = i % internalKmerSize; // gets the starting idx
			}
		}

		for (auto i = orderIdx; i < kmers.size(); i += internalKmerSize)
		{
			optimalKmers.push_back(kmers[i]);
		}

		return optimalKmers;
	}

	std::vector< InternalKmer > KmerLookup::getOptimalKmerSubsetRandom(const std::vector< InternalKmer >& kmers)
	{
		std::vector< InternalKmer > optimalKmers;
		return optimalKmers;
	}
} // namespace dmp
