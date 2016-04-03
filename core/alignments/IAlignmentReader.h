#ifndef DMP_ALIGNMENTS_IALIGNMENTREADER_HPP
#define DMP_ALIGNMENTS_IALIGNMENTREADER_HPP

#include "containers/KmerSetManager.hpp"

#include <memory>

#include <boost/noncopyable.hpp>

namespace dmp
{
	class IAlignmentReader : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentReader > SharedPtr;
		IAlignmentReader() {}
		virtual ~IAlignmentReader() {}

		virtual void processAllReadsInRegion(KmerSetManager::SharedPtr kmerSetManager) = 0;
	};
}

#endif //DMP_IALIGNMENTREADER_HPP
