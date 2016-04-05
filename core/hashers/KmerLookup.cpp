#include "KmerLookup.h"

#include <iostream>

namespace dmp
{
	KmerLookup* KmerLookup::s_instance = NULL;
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
} // namespace dmp
