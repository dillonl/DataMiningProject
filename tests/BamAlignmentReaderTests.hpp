#ifndef TESTS_BAMALIGNMENTREADERTESTS_HPP
#define TESTS_BAMALIGNMENTREADERTESTS_HPP

#include "config/DataConfig.h"
#include "alignments/BamAlignmentReader.h"
#include "alignments/AlignmentRegistration.h"

namespace
{
	using namespace dmp;

	TEST(BamAlignmentReaderTest, GetAllRegions)
	{
		std::string path = NA12878_BAM;
		auto bamAlignmentReader = BamAlignmentReader::CreateSharedPtr(path);
		bamAlignmentReader->processAllReadsInBam();
		std::cout << "done reading" << std::endl;
		AlignmentRegistration::Instance()->AggregateAlignmentsThreads(50);
		std::cout << "count: " << AlignmentRegistration::Instance()->getSize() << std::endl;
	}
}


#endif //TESTS_BAMALIGNMENTREADERTESTS_HPP
