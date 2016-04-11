#ifndef DMP_ALIGNMENTS_ALIGNMENT_H
#define DMP_ALIGNMENTS_ALIGNMENT_H

#include <memory>
#include <unordered_map>
#include <mutex>
#include <vector>
#include <string>

#include "IAlignment.h"
#include "utils/Types.h"

#include <boost/noncopyable.hpp>

namespace dmp
{

	class Alignment : public IAlignment
	{
	public:
		typedef std::shared_ptr< Alignment > SharedPtr;
		static Alignment::SharedPtr CreateAlignment(position pos, bool aligned, const std::vector< InternalKmer >& optimalKmers)
		{
			return std::make_shared< Alignment >(pos, aligned, optimalKmers);
		}

		position getPosition() override { return this->m_position; }
		std::vector< InternalKmer > getOptimalKmerSubsets() override { return this->m_optimal_kmer_subsets; }

	    Alignment(position pos, bool aligned, const std::vector< InternalKmer > optimalKmers) :
            m_position(pos),
			m_is_aligned(aligned),
	        m_optimal_kmer_subsets(optimalKmers)
		{
		}
        ~Alignment() {}

	private:
		std::vector< InternalKmer > m_optimal_kmer_subsets;
		position m_position;
		bool m_is_aligned;
	};
}

#endif //DMP_ALIGNMENTS_ALIGNMENT_H
