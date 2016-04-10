#include "AlignmentRegistration.h"

#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <vector>

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

	void AlignmentRegistration::AggregateAlignmentsThreads(double matchPercent)
	{
		ThreadPool tp;

		std::deque< std::shared_ptr< std::future< void > > > futureFunctions;
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
				continue;
			}
			else
			{
				futureFunctions.emplace_back(futureFunct);
			}
		}
	}

	void AlignmentRegistration::AggregateAlignments(double matchPercent, InternalKmer kmer)
	{
		size_t counter = 0;
		std::unordered_set< InternalKmer > alignmentKmerSet;
		std::vector< IAlignment::SharedPtr > alignmentList = this->m_alignment_ptrs[kmer];
		for (size_t j = 0; j < alignmentList.size(); ++j)
		{
			auto a1SubSet = alignmentList[j]->getOptimalKmerSubsets();
			for (size_t k = (j + 1); k < alignmentList.size(); ++k)
			{
				auto a2SubSet = alignmentList[k]->getOptimalKmerSubsets();
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
				std::unordered_set< InternalKmer > uSet(a1SubSet.begin(), a1SubSet.end());
				uSet.insert(a2SubSet.begin(), a2SubSet.end());

				double percent = ((double)longestRun / (double) uSet.size());
				counter = (percent >= matchPercent) ? counter + 1 : counter;
				if (percent > 20)
				{
					std::cout << "p: " << percent << " c: " << counter << " l: " << longestRun << std::endl;
				}
			}
		}
		std::cout << "processed: " << kmer << std::endl;
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
}// namespace dmp
