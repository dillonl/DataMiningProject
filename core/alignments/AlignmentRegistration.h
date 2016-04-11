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

		uint32_t AggregateAlignmentsThreads(double matchPercent);
		void RegisterAlignment(IAlignment::SharedPtr alignmentPtr);
		uint32_t AggregateAlignments(double matchPercent, InternalKmer kmer);
		size_t getSize();
		std::vector< std::tuple< IAlignment::SharedPtr, IAlignment::SharedPtr, double > > getAlignmentsAndPercent();

	private:
		static AlignmentRegistration* s_instance;
		std::mutex m_lock;
		std::vector< std::vector< IAlignment::SharedPtr > > m_alignment_ptrs;
		std::vector< std::tuple< IAlignment::SharedPtr, IAlignment::SharedPtr, double > > m_alignments_percent;
		std::vector< uint32_t > m_alignment_distances;
	};
} // namespace dmp

#endif //CORE_ALIGNMENT_ALIGNMENTREGISTRATION_H
