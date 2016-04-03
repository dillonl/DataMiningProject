#ifndef DMP_CORE_PARAMETERS_H
#define DMP_CORE_PARAMETERS_H

#include <memory>

#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>

namespace dmp
{
	class Parameters : private boost::noncopyable
	{
	public:
		Parameters();
		~Parameters();

		void parseDMP(int argc, char** argv);
		bool showHelp();
		void printHelp();
		bool validateRequired();

	private:
		std::shared_ptr< boost::program_options::options_description > m_options_description_ptr;
		boost::program_options::variables_map m_variables_map;
	};
}

#endif //DMP_CORE_PARAMETERS_H
