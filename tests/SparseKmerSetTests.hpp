#ifndef TESTS_SPARSEKMERSETTESTS_H
#define TESTS_SPARSEKMERSETTESTS_H

#include <string>
#include <unordered_set>

// #include "containers/KmerSet.hpp"
#include "containers/SparseKmerSet.hpp"
#include "containers/KmerSetManager.hpp"
#include "parsers/AlignmentParser.hpp"
#include "utils/ThreadPool.hpp"

namespace
{
	using namespace rufus;

	TEST(SparseKmerSetTest, InsertionOfZeroItemsSizeZero)
	{
		SparseKmerSet::SharedPtr sks = std::make_shared< SparseKmerSet >();
		ASSERT_EQ(sks->getSetSize(), 0);
	}

	TEST(SparseKmerSetTest, KmerCountZero)
	{
		SparseKmerSet::SharedPtr sks = std::make_shared< SparseKmerSet >();;
		InternalKmer ik = 0; // represents a 25-mer of all A's
		ASSERT_EQ(sks->getKmerCount(ik), 0);
	}

	TEST(SparseKmerSetTest, KmerCountOne)
	{
		SparseKmerSet::SharedPtr sks = std::make_shared< SparseKmerSet >();;
		InternalKmer ik = 0; // represents a 25-mer of all A's
		sks->addKmer(ik);
		ASSERT_EQ(sks->getKmerCount(ik), 1);
	}

	/*
    TEST(AlignmentParserTest, TestInsertionSpeed)
	{
		bool success = true;
		SparseKmerSet kSet;
		uint64_t counter = 0;
		for (uint64_t i = 0; i < 10; ++i)
		{
			for (uint64_t j = 0; j < 5000000; ++j)
			{
				kSet.addKmer((InternalKmer)j);
				++counter;
			}
		}
		std::cout << "total: " << counter << std::endl;
		std::cout << "kmer set size: " << kSet.getSetSize() << std::endl;
		std::cout << "kmer 0 count: " << kSet.getKmerCount(0) << std::endl;
		ASSERT_TRUE(success); //(34588 ms total)
	}
	*/

/*
	TEST(AlignmentParserTest, TestInsertionSpeed2)
	{
		bool success = true;
		KmerSetManager kSetManager;
		static std::atomic< uint64_t > counter(0);
		// auto funct = std::bind(&KmerSetManager::addKmer, &kSetManager, std::placeholders::_1);
		auto funct = [&] ()
		{
			for (uint64_t i = 0; i < 2500000; ++i)
			{
				kSetManager.addKmer((InternalKmer)counter++);
			}
			std::cout << "finished one loop" << std::endl;
		};
		for (int i = 0; i < 5; ++i)
		{
			ThreadPool::Instance()->enqueue(funct);
			ThreadPool::Instance()->enqueue(funct);
			ThreadPool::Instance()->enqueue(funct);
			ThreadPool::Instance()->enqueue(funct);
		}
		kSetManager.stop();
		ThreadPool::Instance()->joinAll();
		auto kSet = kSetManager.getKmerSet();
		std::cout << "kmer count: " << kSet->getSetSize() << std::endl;
		std::cout << "total count: " << counter << std::endl;
		ASSERT_TRUE(success); //(72470 ms total)
	}
*/

	/*
	TEST(AlignmentParserTest, ParseA)
	{
		//std::string alignment = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		std::string alignment = "CCCCCCCCCCCCCCCCCCCCCCCCCA"; //AAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTT";
		InternalKmer parsedKmer;

		std::unordered_set< InternalKmer, KmerHash, KmerKeyEqual > kSet;
		ThreadPool::Instance()->setThreadCount(5);
		std::cout << "starting" << std::endl;
		bool success = true;
		for (uint32_t i = 0; i < 1; ++i)
		{
			ThreadPool::Instance()->enqueue(std::bind(&AlignmentParser::ParseAlignment< 25 >, alignment.c_str(), alignment.size(), kSet));
			//success |= AlignmentParser::ParseAlignment(alignment.c_str(), alignment.size(), parsedKmer);
		}
		ThreadPool::Instance()->joinAll();
		std::cout << "ending" << std::endl;
		ASSERT_TRUE(success);
		//ASSERT_EQ();
	}
	*/

	/*
	TEST(AlignmentParserTest, HashTest)
	{
		//std::string alignment = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		// std::string alignment = "CCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTT";
		std::string alignment = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
		InternalKmer parsedKmer;

		std::unordered_set< InternalKmer, KmerHash, KmerKeyEqual > kSet;
		ThreadPool::Instance()->setThreadCount(1);
		std::cout << "starting" << std::endl;
		bool success = true;
		for (uint32_t i = 0; i < 2; ++i)
		{
			ThreadPool::Instance()->enqueue(std::bind(&AlignmentParser::ParseAlignment, alignment.c_str(), alignment.size(), &kSet));
			//success |= AlignmentParser::ParseAlignment(alignment.c_str(), alignment.size(), parsedKmer);
		}
		ThreadPool::Instance()->joinAll();
		std::cout << "ending: " << kSet.size() << std::endl;
		ASSERT_TRUE(success);
		//ASSERT_EQ();
	}
	*/
}

#endif //TESTS_SPARSEKMERSETTESTS_H
