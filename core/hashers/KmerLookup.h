#ifndef DMP_KMERLOOKUP_HPP
#define DMP_KMERLOOKUP_HPP

#include "utils/Types.h"

#include <vector>

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
			}
			return s_instance;
		}

		std::vector< InternalKmer > getOptimalKmerSubset(const std::vector< InternalKmer >& kmers);
	private:
		KmerLookup() {}
		~KmerLookup() {}

		static KmerLookup* s_instance;
	};
} // namespace dmp

#endif // DMP_KMERLOOKUP_HPP
