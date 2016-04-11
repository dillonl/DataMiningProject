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

		static inline std::string UnParseAlignment(std::vector< InternalKmer >& kmers)
		{
			static std::unordered_map< InternalKmer, char > kmerMap = {{0,'A'}, {1,'C'}, {2, 'T'}, {3, 'G'}};
			std::string alignment;
			for (auto& kmer : kmers)
			{
				for (uint32_t i = 0; i < sizeof(kmer); ++i)
				{
					uint32_t shifter = i*8;
					auto k1 = 0x3 & (kmer >> shifter);
					auto k2 = 0x3 & (kmer >> shifter + 2);
					auto k3 = 0x3 & (kmer >> shifter + 4);
					auto k4 = 0x3 & (kmer >> shifter + 6);
					alignment += kmerMap[k1];
					alignment += kmerMap[k2];
					alignment += kmerMap[k3];
					alignment += kmerMap[k4];
				}
			}
			return alignment;
		}

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
