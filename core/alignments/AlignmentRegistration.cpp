#include "AlignmentRegistration.h"

#include <unordered_set>
#include <iostream>

namespace dmp
{
	AlignmentRegistration* AlignmentRegistration::s_instance = NULL;

	AlignmentRegistration::AlignmentRegistration()
	{
		this->m_alignment_ptrs.resize((size_t)std::numeric_limits< InternalKmer >::max());
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
			// if (!std::get< 1 >(alignmentKmerSet.emplace(kmer))) { continue; }
			this->m_alignment_ptrs[(uint32_t)kmer].emplace_back(alignmentPtr);
		}
	}

	void AlignmentRegistration::AggregateAlignments(double matchPercent)
	{
		std::lock_guard< std::mutex > lock(this->m_lock);
		std::unordered_set< InternalKmer > alignmentKmerSet;
		for (uint32_t i = 0; i < std::numeric_limits< InternalKmer >::max(); ++i)
		{
			std::vector< IAlignment::SharedPtr > alignmentList = this->m_alignment_ptrs[i];
		}
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
