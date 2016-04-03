INCLUDE_DIRECTORIES(
		${SCI_BOOST_INCLUDE}
		${BAMTOOLS_INCLUDE}
		)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

set(DMP_UTIL_SOURCES
  utils/Parameters.cpp
  )

set(DMP_ALIGNMENTS_SOURCES
  alignments/BamAlignmentReader.cpp
  )

add_library(dmp_core
  ${DMP_UTIL_SOURCES}
  ${DMP_ALIGNMENTS_SOURCES}
  )

SET(CORE_LIB dmp_core CACHE STRING "THE CORE LIBRARY") #add dmp_core and dmp_utils to the list of core libs


add_dependencies(${CORE_LIB} ${BAMTOOLS_PROJECT} Boost_external )
