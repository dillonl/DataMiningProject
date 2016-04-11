#include "AlignmentRegistration.h"
#include "parsers/AlignmentParser.hpp"

#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <vector>
#include <tuple>

namespace dmp
{
	AlignmentRegistration* AlignmentRegistration::s_instance = NULL;

	AlignmentRegistration::AlignmentRegistration()
	{
		this->m_alignment_ptrs.resize((size_t)std::numeric_limits< InternalKmer >::max() + 1);
	}

	AlignmentRegistration::~AlignmentRegistration()
	{
	}

	void AlignmentRegistration::RegisterAlignment(IAlignment::SharedPtr alignmentPtr)
	{
		std::lock_guard< std::mutex > lock(this->m_lock);
		std::unordered_set< InternalKmer > alignmentKmerSet;
		for (auto kmer : alignmentPtr->getOptimalKmerSubsets())
		{
			if (!std::get< 1 >(alignmentKmerSet.emplace(kmer))) { continue; }
			this->m_alignment_ptrs[(uint32_t)kmer].emplace_back(alignmentPtr);
		}
	}

	uint32_t AlignmentRegistration::AggregateAlignmentsThreads(double matchPercent)
	{
		uint32_t totalCount = 0;
		ThreadPool tp;

		std::deque< std::shared_ptr< std::future< uint32_t > > > futureFunctions;
		size_t counter = 0;
		std::unordered_set< InternalKmer > alignmentKmerSet;
		for (InternalKmer i = 0; i < std::numeric_limits< InternalKmer >::max(); ++i)
		{
			auto funct = std::bind(&AlignmentRegistration::AggregateAlignments, this, matchPercent, i);
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

	uint32_t AlignmentRegistration::AggregateAlignments(double matchPercent, InternalKmer kmer)
	{
		// std::lock_guard< std::mutex > l(m_lock);
		size_t counter = 0;
		std::unordered_set< InternalKmer > alignmentKmerSet;
		std::vector< IAlignment::SharedPtr > alignmentList = this->m_alignment_ptrs[kmer];
		for (size_t j = 0; j < alignmentList.size(); ++j)
		{
			auto a1 = alignmentList[j];
			auto a1SubSet = a1->getOptimalKmerSubsets();
			for (size_t k = (j + 1); k < alignmentList.size(); ++k)
			{
				auto a2 = alignmentList[k];
				auto a2SubSet = a2->getOptimalKmerSubsets();

				std::vector< InternalKmer > intersection(a1SubSet.size());
				std::vector< InternalKmer > unionSet(a1SubSet.size() + a2SubSet.size());
				std::unordered_set< InternalKmer > a1Set(a1SubSet.begin(), a1SubSet.end());
				std::unordered_set< InternalKmer > a2Set(a2SubSet.begin(), a2SubSet.end());
				auto iter1 = std::set_intersection(a1Set.begin(), a1Set.end(), a2Set.begin(), a2Set.end(), intersection.begin());
				intersection.resize(iter1-intersection.begin());
				auto iter2 = std::set_union(a1Set.begin(), a1Set.end(), a2Set.begin(), a2Set.end(), unionSet.begin());
				unionSet.resize(iter2-unionSet.begin());

				auto jacobScore = (double)intersection.size() / (double)unionSet.size();
				if (jacobScore < matchPercent) { continue; }

				auto consecutiveRun = 0;
				auto longestRun = 0;
				for (auto x = 0; x < a1SubSet.size(); ++x)
				{
					for (auto y = 0; y < a2SubSet.size(); ++y)
					{
						if (a1SubSet[x] == a2SubSet[y])
						{
							++consecutiveRun;
							++x;
							if (x >= a1SubSet.size()) { break; }
						} else {
							longestRun = (consecutiveRun > longestRun) ? consecutiveRun : longestRun;
							consecutiveRun = 0;
						}
					}
				}
				longestRun = (consecutiveRun > longestRun) ? consecutiveRun : longestRun;

				double percent = ((double)longestRun / (double) a1SubSet.size());
				if (percent >= matchPercent)
				{
					++counter;
					m_alignments_percent.emplace_back(std::make_tuple(a1, a2, percent));
				}
			}
		}
		return counter;
	}

	size_t AlignmentRegistration::getSize()
	{
		std::lock_guard< std::mutex > lock(this->m_lock);
		std::unordered_set< InternalKmer > alignmentKmerSet;
		size_t size = 0;
		for (uint32_t i = 0; i < std::numeric_limits< InternalKmer >::max(); ++i)
		{
			size += this->m_alignment_ptrs[i].size();
		}
		return size;
	}

	std::vector< std::tuple< IAlignment::SharedPtr, IAlignment::SharedPtr, double > > AlignmentRegistration::getAlignmentsAndPercent()
	{
		return this->m_alignments_percent;
	}
}// namespace dmp
