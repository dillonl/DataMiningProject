#ifndef DMP_RESULTSPRINTER_HPP
#define DMP_RESULTSPRINTER_HPP

#include "alignments/AlignmentRegistration.h"

#include <iostream>
#include <fstream>

namespace dmp
{
	class ResultsPrinter
	{
	public:
		static void printResults(const std::string& outPrefix, AlignmentRegistration* alignmentRegistration)
		{
			std::ofstream distanceFile;
			std::ofstream percentFile;

			distanceFile.open(outPrefix + ".distance.txt");
			percentFile.open(outPrefix + ".percent.txt");

			size_t counter = 0;
			size_t elementSize = alignmentRegistration->getAlignmentsAndPercent().size();
			for (auto alignmentsAndPercent : alignmentRegistration->getAlignmentsAndPercent())
			{
				IAlignment::SharedPtr a1 = std::get< 0 >(alignmentsAndPercent);
				IAlignment::SharedPtr a2 = std::get< 1 >(alignmentsAndPercent);
				auto percent = std::get< 2 >(alignmentsAndPercent);

				uint32_t dist = (a1->getPosition() < a2->getPosition()) ? a2->getPosition() - a1->getPosition() : a1->getPosition() - a2->getPosition();
				distanceFile << dist << std::endl;
				percentFile << percent << std::endl;
				// std::cout << "dist: " << dist << " perc: " << percent << " p: " << a1->getPosition() << " " << a2->getPosition() << std::endl;
				// oFile <<
			}
			distanceFile.close();
			percentFile.close();
			// oFile.close();
		}
	};
}

#endif //DMP_RESULTSPRINTER_HPP
