#ifndef CORE_PARSERS_ALIGNMENTPARSER_H
#define CORE_PARSERS_ALIGNMENTPARSER_H

#include "utils/Types.h"

#include <memory>
#include <iostream>
#include <bitset>
#include <unordered_map>

#include <boost/noncopyable.hpp>

namespace dmp
{

	class AlignmentParser : private boost::noncopyable
	{
	public:
		static inline bool ParseAlignment(const char* alignment, size_t kmerIterations, std::vector< InternalKmer >& kmers)
		{
			for (auto i = 0; i < kmerIterations; ++i)
			{
				kmers[i] = 0;
				bool validKmer = true;
				validKmer &= unroller(i, uint_< KMER_SIZE - 1 >(),  alignment, kmers[i]);
			}
			return true;
		}

		/*
		static inline bool ParseAlignment(const char* alignment, size_t kmerIterations, InternalKmer* kmers)
			{
				for (auto i = 0; i < kmerIterations; ++i)
				{
					InternalKmer internalKmer = 0; // we must guarentee internalKmer is zero
					bool validKmer = true;
					validKmer &= unroller(i, uint_< KMER_SIZE - 1 >(),  alignment, internalKmer);
					kmers[i] = internalKmer;
				}
				return true;
			}
		*/

	private:
		static const InternalKmer ShiftByOne = 1; // these need to the same as InternalKmer
		static const InternalKmer ShiftByThree = 3;

		AlignmentParser() = delete;
		~AlignmentParser() = delete;

		template <InternalKmer N> struct uint_{ };

		template <InternalKmer N, typename IterT>
		static inline bool unroller(const IterT& iter, uint_<N>, const char* alignment, InternalKmer& internalKmer)
		{
			unroller(iter, uint_<N-1>(), alignment, internalKmer);
			return parse(alignment, internalKmer, N, iter);
		}

		template <typename IterT>
		static inline bool unroller(const IterT& iter, uint_<0>, const char* alignment, InternalKmer& internalKmer)
		{
			return parse(alignment, internalKmer, 0, iter);
		}

		static inline bool parse(const char* alignmentCharacter, InternalKmer& kmer, InternalKmer kmerBitOffset, uint64_t counter)
		{
			InternalKmer shifter = (kmerBitOffset * 2);
			if ((alignmentCharacter[kmerBitOffset + counter]) & (1 >> (3))) { return false; } // if this is a non basepair char
			kmer |= (((alignmentCharacter[kmerBitOffset + counter] >> ShiftByOne) & 0x3) & ShiftByThree) << shifter;
			return true;
		}
	};
}

#endif //CORE_PARSERS_ALIGNMENTPARSER_H
