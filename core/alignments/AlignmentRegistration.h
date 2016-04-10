#ifndef CORE_ALIGNMENT_ALIGNMENTREGISTRATION_H
#define CORE_ALIGNMENT_ALIGNMENTREGISTRATION_H

#include "IAlignment.h"
#include "utils/ThreadPool.hpp"

#include <boost/noncopyable.hpp>

namespace dmp
{
	class AlignmentRegistration : boost::noncopyable
	{
	public:
		AlignmentRegistration();
		~AlignmentRegistration();

		static AlignmentRegistration* Instance()
		{
			if (s_instance == NULL)
			{
				s_instance = new AlignmentRegistration();
			}
			return s_instance;
		}

		void AggregateAlignmentsThreads(double matchPercent);
		void RegisterAlignment(IAlignment::SharedPtr alignmentPtr);
		void AggregateAlignments(double matchPercent, InternalKmer kmer);
		size_t getSize();
	private:
		static AlignmentRegistration* s_instance;
		std::mutex m_lock;
		std::vector< std::vector< IAlignment::SharedPtr > > m_alignment_ptrs;
	};
} // namespace dmp

#endif //CORE_ALIGNMENT_ALIGNMENTREGISTRATION_H
