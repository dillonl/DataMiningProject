#ifndef DMP_ALIGNMENTS_IALIGNMENT_H
#define DMP_ALIGNMENTS_IALIGNMENT_H

#include <memory>
#include <unordered_map>
#include <mutex>
#include <vector>
#include <string>

#include "utils/Types.h"

#include <boost/noncopyable.hpp>

namespace dmp
{

	class IAlignment : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignment > SharedPtr;

		virtual position getPosition() = 0;
		virtual std::vector< InternalKmer > getOptimalKmerSubsets() = 0;

	};
}

#endif //DMP_ALIGNMENTS_IALIGNMENT_H
