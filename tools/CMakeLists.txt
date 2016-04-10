# =================================
# dmp
#
# tools/CMakeLists.txt
# =================================

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

# set header and source files
set(DMP_TOOLS_SOURCES
	dmp.cpp
)

INCLUDE_DIRECTORIES(
  ${SCI_BOOST_INCLUDE}
  ${BAMTOOLS_INCLUDE}
  )


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # clang Doesnt use pthread
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

#set(CMAKE_BUILD_TYPE Debug) # uncomment this when you need a debug build

#create binary
add_executable(dmp
  ${DMP_TOOLS_SOURCES}
)

add_dependencies(dmp ${BAMTOOLS_PROJECT} Boost_external)

# link libraries to executables
target_link_libraries(dmp
    ${BAMTOOLS_LIB}
   	${BAMTOOLS_UTIL_LIB}
	${ZLIB_LIB}
	${CORE_LIB}
	#${FASTAHACK_LIB}
	${SCI_BOOST_LIBRARY}
)

