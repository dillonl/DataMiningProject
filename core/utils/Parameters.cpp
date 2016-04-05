#include "Parameters.h"

#include <thread>
#include <iostream>

namespace dmp
{
	Parameters::Parameters()
	{
	}

	Parameters::~Parameters()
	{
	}

	void Parameters::parseDMP(int argc, char** argv)
	{
		m_options_description_ptr = std::make_shared< boost::program_options::options_description >("options");
		m_options_description_ptr->add_options()
			("help,h","Print help message")
			(",b", boost::program_options::value< std::string >()->default_value(NULL), "Path to input BAM file")
			(",o", boost::program_options::value< std::string >()->default_value(""), "Path to output directory [optional - if not provided then prints to std::cout]")
			(",t", boost::program_options::value< uint32_t >()->default_value(std::thread::hardware_concurrency() * 2), "Thread count [optional - default is number of cores x 2]");
		auto parseCommandLine = boost::program_options::parse_command_line(argc, argv, *m_options_description_ptr);
		boost::program_options::store(parseCommandLine, m_variables_map);
	}

	bool Parameters::showHelp()
	{
		return (m_variables_map.count("help") > 0);
	}

	void Parameters::printHelp()
	{
		std::cout << *m_options_description_ptr << std::endl;
	}

	bool Parameters::validateRequired()
	{
		try
		{
			boost::program_options::notify(m_variables_map);
		}
		catch (boost::program_options::required_option& e)
		{
			std::cout << "ERROR " << e.what() << std::endl << std::endl;
			return false;
		}
		return true;
	}
}
