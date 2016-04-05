#include "gtest/gtest.h"

//#include "config/DataPaths.h"
#include "BamAlignmentReaderTests.hpp"
//#include "AlignmentParserTests.hpp"
//#include "SparseKmerSetTests.hpp"

GTEST_API_ int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
