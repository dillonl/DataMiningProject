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


# Setting up external library for BWA
SET(BWA_PROJECT bwa_project CACHE INTERNAL "bwa project name")
SET(BWA_DIR ${CMAKE_BINARY_DIR}/externals/bwa CACHE INTERNAL "bwa project directory")
ExternalProject_Add(${BWA_PROJECT}
	GIT_REPOSITORY https://github.com/lh3/bwa.git
	GIT_TAG 0911122749ccfa7bd4ae9c464778f56372524f96 #lock in the commit id so we don't this doesn't break in the future
	UPDATE_COMMAND ""
	INSTALL_COMMAND ""
	BUILD_COMMAND "make"
	CONFIGURE_COMMAND ""
	BINARY_DIR ${BWA_DIR}
	SOURCE_DIR ${BWA_DIR}
	DEPENDS ${ZLIB_PROJECT}
    CMAKE_CACHE_ARGS
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}

)

SET(BWA_LIB ${BWA_DIR}/libbwa.a CACHE INTERNAL "BWA Library")
SET(BWA_INCLUDE ${BWA_DIR} CACHE INTERNAL "BWA Include")
