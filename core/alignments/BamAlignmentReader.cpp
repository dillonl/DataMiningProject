#include "BamAlignmentReader.h"

#include "containers/KmerSet.hpp"
#include "hashers/KmerLookup.h"
#include "Alignment.h"
#include "AlignmentRegistration.h"
#include "parsers/AlignmentParser.hpp"
#include "utils/ThreadPool.hpp"

#include <deque>
#include <limits>
#include <cstdlib>

namespace dmp
{
    BamAlignmentReader::SharedPtr BamAlignmentReader::CreateSharedPtr(const std::string& filePath)
	{
        return std::make_shared< BamAlignmentReader >(filePath);
    }

    BamAlignmentReader::BamAlignmentReader(const std::string& filePath) :
		m_file_path(filePath)
	{
	}

	BamAlignmentReader::~BamAlignmentReader()
	{
	}

	std::vector< BamAlignmentReader::BamRegion::SharedPtr > BamAlignmentReader::getAllSpacedOutRegions()
	{
		BamTools::BamReader bamReader;
		if (!bamReader.Open(this->m_file_path))
		{
			throw "Unable to open bam file";
		}
		bamReader.LocateIndex();

		// get all region ids
		std::vector< int > regionIDs;
		for (auto refData : bamReader.GetReferenceData())
		{
			regionIDs.emplace_back(bamReader.GetReferenceID(refData.RefName));
		}

		auto referenceData = bamReader.GetReferenceData();
		// get the region pointers
		std::vector< BamRegion::SharedPtr > regionPtrs;
		uint32_t intervalSize = 1000000;
		for (auto regionID : regionIDs)
		{
			uint32_t regionLastPosition = referenceData[regionID].RefLength;
			uint32_t currentPosition = 0;
			while (currentPosition < regionLastPosition)
			{
				uint32_t positionDelta = ((currentPosition + intervalSize) > regionLastPosition) ? (regionLastPosition - currentPosition) : intervalSize - 1;
				uint32_t endPosition = currentPosition + positionDelta;
				auto bamRegionPtr = std::make_shared< BamRegion >(regionID, currentPosition, endPosition);
				regionPtrs.emplace_back(bamRegionPtr);
				currentPosition += positionDelta + 1;
				// bamRegionPtr->print();
			}
		}

		bamReader.Close();
		return regionPtrs;
	}

	uint32_t BamAlignmentReader::processAllReadsInBam()
	{
		uint32_t totalCount = 0;
		ThreadPool tp;
		auto spacedOutRegions = getAllSpacedOutRegions();
		std::deque< std::shared_ptr< std::future< uint32_t > > > futureFunctions;
		for (auto regionPtr : spacedOutRegions)
		{
			auto funct = std::bind(&BamAlignmentReader::processReads, this, regionPtr);
			auto futureFunct = tp.enqueue(funct);
			futureFunctions.emplace_back(futureFunct);
		}
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready)
			{
				totalCount += futureFunct->get();
				continue;
			}
			else
			{
				futureFunctions.emplace_back(futureFunct);
			}
		}
		return totalCount;
	}

	uint32_t BamAlignmentReader::processReads(BamRegion::SharedPtr bamRegionPtr)
	{
		uint32_t counter = 0;
		BamTools::BamReader bamReader;
		if (!bamReader.Open(this->m_file_path))
		{
			throw "Unable to open bam file";
		}
		bamReader.LocateIndex();
		bamReader.SetRegion(bamRegionPtr->getRegionID(), bamRegionPtr->getStartPosition(), bamRegionPtr->getRegionID(), bamRegionPtr->getEndPosition());

		auto bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		std::vector< InternalKmer > internalKmers(100);
		size_t kmerCount = 0;

		uint32_t count = 0;

		while(bamReader.GetNextAlignment(*bamAlignmentPtr))
		{
			if (bamAlignmentPtr->Position < bamRegionPtr->getStartPosition()) { continue; }
			auto kmersNumber = (bamAlignmentPtr->Length - KMER_SIZE);
			if (kmersNumber > internalKmers.size()) { internalKmers.resize(kmersNumber); }
			if (AlignmentParser::ParseAlignment(bamAlignmentPtr->QueryBases.c_str(), kmersNumber, internalKmers))
			{
				auto optimalKmerSubsets = KmerLookup::Instance()->getOptimalKmerSubset(internalKmers);
				// auto optimalKmerSubsets = KmerLookup::Instance()->getOptimalKmerSubsetOrdered(internalKmers);
				auto alignmentPtr = std::make_shared< Alignment >(bamAlignmentPtr->Position, bamAlignmentPtr->IsMapped(), optimalKmerSubsets);
				// auto alignmentPtr = Alignment::CreateAlignment(bamAlignmentPtr->Position, );
				AlignmentRegistration::Instance()->RegisterAlignment(alignmentPtr);
			}
			++count;
			// if (count >= 50) { break; }
		}
		bamReader.Close();
		return count;

		// static std::mutex l;
		// l.lock();
		// std::cout << "counter: " << sCounter << std::endl;
		// l.unlock();
		// std::cout << "total count: [" << m_kmer_set_ptr->getSetSize() << "] " << counter << " ";
        // bamRegionPtr->print();
	}
}
