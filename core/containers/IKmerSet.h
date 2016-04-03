#ifndef CORE_CONTAINERS_IKMERSET_HPP
#define CORE_CONTAINERS_IKMERSET_HPP

#include "utils/Types.h"

#include <memory>

#include <boost/noncopyable.hpp>

namespace dmp
{
	class IKmerSet : boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IKmerSet > SharedPtr;
		IKmerSet()
		{
		}

		virtual ~IKmerSet()
		{
		}


		virtual void resize(uint64_t size) = 0;
		virtual void addKmer(InternalKmer internalKmer) = 0;
		virtual uint64_t getKmerCount(InternalKmer internalKmer) = 0;
		virtual size_t getSetSize() = 0;
		virtual void addAllKmersToPassedInSet(IKmerSet::SharedPtr kmerSetPtr) = 0;
	};
}

#endif // CORE_CONTAINERS_IKMERSET_HPP
