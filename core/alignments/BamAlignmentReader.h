#ifndef DMP_ALIGNMENTS_BAMALIGNMENTREADER_HPP
#define DMP_ALIGNMENTS_BAMALIGNMENTREADER_HPP

#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "IAlignmentReader.h"
//#include "containers/KmerSetManager.hpp"

#include <boost/noncopyable.hpp>

#include <unordered_set>

namespace dmp
{
	class BamAlignmentReader : private boost::noncopyable
	{
	private:
		class BamRegion : private boost::noncopyable
	    {
		public:
			typedef std::shared_ptr< BamRegion > SharedPtr;
		BamRegion(int regionID, int startPosition, int endPosition) :
			m_region_id(regionID),
				m_start_position(startPosition),
				m_end_position(endPosition)
				{
				}

			int getRegionID() { return this->m_region_id; }
			int getStartPosition() { return this->m_start_position; }
			int getEndPosition() { return this->m_end_position; }
			void print() { std::cout << m_region_id << " " << m_start_position << " " << m_end_position << std::endl; }
		private:
			int m_region_id;
			int m_start_position;
			int m_end_position;
		};

	public:
        typedef std::shared_ptr< BamAlignmentReader > SharedPtr;
		static BamAlignmentReader::SharedPtr CreateSharedPtr(const std::string& filePath);
        BamAlignmentReader(const std::string& filePath);
		~BamAlignmentReader();

		uint32_t processAllReadsInBam();
		/* void processAllReadsInRegion(SparseKmerSet::SharedPtr kmerSetPtr); */

		/* static std::vector< int > getAllRegionsInBam(const std::string& filePath); */

	private:
		uint32_t processReads(BamRegion::SharedPtr bamRegionPtr);
		std::vector< BamRegion::SharedPtr > getAllSpacedOutRegions();

		std::string m_file_path;
		std::mutex m_lock;
		/* std::unordered_set< InternalKmer, KmerHash, KmerKeyEqual > m_set; */
	};
}

#endif //DMP_BAMALIGNMENTREADER_HPP
