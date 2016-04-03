#include "BamAlignmentReader.h"

#include "containers/KmerSet.hpp"
#include "parsers/AlignmentParser.hpp"
#include "utils/ThreadPool.hpp"

#include <deque>
#include <limits>
#include <cstdlib>

namespace dmp
{
	BamAlignmentReader::BamAlignmentReader(const std::string& filePath) :
		m_file_path(filePath)
	{
		// m_kmer_mutexes.resize(256);
		for (uint32_t i = 0; i < 256; ++i)
		{
			m_kmer_set_ptrs.emplace_back(std::make_shared< KmerSet >());
			m_kmer_mutexes.emplace_back(std::make_shared< std::mutex >());
		}
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
				// if (endPosition != regionLastPosition && endPosition + 10000 > regionLastPosition) { endPosition = regionLastPosition; } // so we don't end up with small regions near the end of the bamregion
				auto bamRegionPtr = std::make_shared< BamRegion >(regionID, currentPosition, endPosition);
				regionPtrs.emplace_back(bamRegionPtr);
				currentPosition += positionDelta + 1;
				// bamRegionPtr->print();
			}
		}

		bamReader.Close();
		return regionPtrs;
	}

	void BamAlignmentReader::processAllReadsInBam()
	{
		ThreadPool tp;
		auto spacedOutRegions = getAllSpacedOutRegions();
		// std::vector< std::shared_ptr< std::future< IKmerSet::SharedPtr > > > futureFunctions;
		std::deque< std::shared_ptr< std::future< IKmerSet::SharedPtr > > > futureFunctions;
		for (auto regionPtr : spacedOutRegions)
		{
			// regionPtr->print();
			auto funct = std::bind(&BamAlignmentReader::processReads, this, regionPtr);
			auto futureFunct = tp.enqueue(funct);
			futureFunctions.emplace_back(futureFunct);
			static int counter = 0;
			// if (counter++ >= 75) { break; }
		}
		// for (auto& futureFunct : futureFunctions)
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready)
			{
				// auto kmerSetPtr = futureFunct->get();
				// std::lock_guard< std::mutex > guard(m_lock);
				// m_set.insert(set.begin(), set.end());
 				// m_kmer_set_ptr->addAllKmersToPassedInSet(kmerSetPtr);
				continue;
			}
			else
			{
				futureFunctions.emplace_back(futureFunct);
			}
			// futureFunct->wait();
		}
	}

	IKmerSet::SharedPtr BamAlignmentReader::processReads(BamRegion::SharedPtr bamRegionPtr)
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

		while(bamReader.GetNextAlignment(*bamAlignmentPtr))
		{
			if (bamAlignmentPtr->Position < bamRegionPtr->getStartPosition()) { continue; }
			auto kmersNumber = (bamAlignmentPtr->Length - KMER_SIZE);
			if (kmersNumber > internalKmers.size()) { internalKmers.resize(kmersNumber); }
			if (AlignmentParser::ParseAlignment(bamAlignmentPtr->QueryBases.c_str(), kmersNumber, internalKmers))
			{
				for (auto i = 0; i < kmersNumber; ++i)
				{
					int idx = internalKmers[i] & 0xff;
					m_kmer_mutexes[idx]->lock();
					m_kmer_set_ptrs[idx]->addKmer(internalKmers[i]);
					m_kmer_mutexes[idx]->unlock();
                    ++counter;
				}
			}
		}
		bamReader.Close();

		// std::cout << "total count: [" << m_kmer_set_ptr->getSetSize() << "] " << counter << " ";
        bamRegionPtr->print();

		// return kmerSetPtr;
		return nullptr;
	}
}
