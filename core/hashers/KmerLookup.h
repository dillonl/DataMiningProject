#ifndef DMP_KMERLOOKUP_HPP
#define DMP_KMERLOOKUP_HPP

#include "utils/Types.h"

#include <vector>
#include <unordered_map>
#include <mutex>

#include <boost/noncopyable.hpp>

namespace dmp
{
	class KmerLookup : private boost::noncopyable
	{
	public:
		static KmerLookup* Instance()
		{
			if (s_instance == NULL)
			{
				s_instance = new KmerLookup();
				s_instance->init();
			}
			return s_instance;
		}

		std::vector< InternalKmer > getOptimalKmerSubset(const std::vector< InternalKmer >& kmers);
		std::vector< InternalKmer > getOptimalKmerSubsetOrdered(const std::vector< InternalKmer >& kmers);
		std::vector< InternalKmer > getOptimalKmerSubsetRandom(const std::vector< InternalKmer >& kmers);
	private:
		KmerLookup();
		~KmerLookup();

		void init();

		static KmerLookup* s_instance;
		std::mutex m_map_init_mutex;
		std::unordered_map< InternalKmer, uint32_t > m_order_mapping;
	};
} // namespace dmp

#endif // DMP_KMERLOOKUP_HPP
