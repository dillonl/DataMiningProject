#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2015 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software. 
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

# Setting up external library openmp

SET(SPARSEHASH_PROJECT sparsehash_project CACHE INTERNAL "sparsehash project name")
SET(SPARSEHASH_DIR ${CMAKE_BINARY_DIR}/externals/sparsehash CACHE INTERNAL "sparsehash project directory")
SET(SPARSEHASH_LIB)
ExternalProject_Add(${SPARSEHASH_PROJECT}
	GIT_REPOSITORY https://github.com/sparsehash/sparsehash.git
	GIT_TAG 804c039c10af1c7a130422735ef7cbf4a7bf2663
    UPDATE_COMMAND ""
    INSTALL_COMMAND "make"
    BUILD_COMMAND ./configure CC=${CMAKE_CXX_COMPILER} CXXFLAGS=-std=c++11
    CONFIGURE_COMMAND ./autogen.sh
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${SPARSEHASH_DIR}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

ExternalProject_Get_Property(${SPARSEHASH_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${SPARSEHASH_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${SPARSEHASH_PROJECT} BINARY_DIR)

#SET(SPARSEHASH_LIB ${SOURCE_DIR}/sparsehash CACHE INTERNAL "sparsehash Lib")
SET(SPARSEHASH_INCLUDE_DIR ${SOURCE_DIR}/src CACHE INTERNAL "sparsehash Include")



