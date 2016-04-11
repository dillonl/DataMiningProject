#ifndef TESTS_BAMALIGNMENTREADERTESTS_HPP
#define TESTS_BAMALIGNMENTREADERTESTS_HPP

#include "config/DataConfig.h"
#include "alignments/BamAlignmentReader.h"
#include "alignments/AlignmentRegistration.h"
#include "utils/ResultsPrinter.hpp"

namespace
{
	using namespace dmp;

	TEST(BamAlignmentReaderTest, GetAllRegions)
	{
		std::string path = NA12878_BAM;
		auto bamAlignmentReader = BamAlignmentReader::CreateSharedPtr(path);
		auto bamReadCount = bamAlignmentReader->processAllReadsInBam();
		std::cout << "done reading" << std::endl;
		auto alignmentCount = AlignmentRegistration::Instance()->AggregateAlignmentsThreads(0.3);
		std::cout << "ac: " << alignmentCount << std::endl;
		std::cout << "rc: " << bamReadCount << std::endl;
		ResultsPrinter::printResults("match_0", AlignmentRegistration::Instance());
		// std::cout << "count: " << AlignmentRegistration::Instance()->getSize() << std::endl;
	}
}


#endif //TESTS_BAMALIGNMENTREADERTESTS_HPP
