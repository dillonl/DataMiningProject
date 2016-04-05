# =================================
# dmp - test
#
# text/CMakeLists.txt
# =================================

# set header and source files
set(DMP_TEST_SOURCES
  gtest_main.cpp
)

# Where Google Test's .h files can be found.
include_directories(
		${CMAKE_SOURCE_DIR}
		${SCI_BOOST_INCLUDE}
		${BAMTOOLS_INCLUDE}
		${GTEST_INCLUDE}
		${GTEST_INCLUDE}/include
  )

set(Boost_USE_MULTITHREADED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()
#set(CMAKE_BUILD_TYPE Debug)


#create binary
add_executable(dmp_tests
  ${DMP_TEST_SOURCES}
  )

add_dependencies(dmp_tests ${BAMTOOLS_PROJECT} Boost_external ${GTEST_PROJECT})


target_link_libraries(dmp_tests
	${BAMTOOLS_LIB}
   	${BAMTOOLS_UTIL_LIB}
	${ZLIB_LIB}
	${CORE_LIB}
	#${FASTAHACK_LIB}
	${SCI_BOOST_LIBRARY}
  	${GTEST_LIB}
  )

