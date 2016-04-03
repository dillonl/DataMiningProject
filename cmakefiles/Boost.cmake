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

# Borrowed CMake code from the MaidSafe Boost CMake build
# found at https://github.com/maidsafe/MaidSafe/blob/master/cmake_modules/add_boost.cmake
# and code borrowed from ITK4 HDFMacros.cmake

SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

# disable auto linking
# also set in Seg3D?
SET(boost_CXX_Flags "-DBOOST_ALL_NO_LIB=1")
IF(APPLE)
  LIST(APPEND boost_CXX_Flag "-DBOOST_LCAST_NO_WCHAR_T" "-DBOOST_THREAD_DONT_USE_ATOMIC")
ENDIF()
IF(WIN32)
  LIST(APPEND boost_CXX_Flag "-DBOOST_BIND_ENABLE_STDCALL")
ENDIF()

# explicitly set library list
SET(boost_Libraries
  "atomic"
  "date_time"
  "exception"
  "filesystem"
  "program_options"
  "regex"
  "serialization"
  "system"
  "thread"
  CACHE INTERNAL "Boost library name.")

# for travis clang builds--need a narrower test
IF(UNIX)
  ADD_DEFINITIONS(-DBOOST_NO_CXX11_ALLOCATOR)
ENDIF()

# TODO: set up 64-bit build detection
# Boost Jam needs to have 64-bit build explicitly configured
IF(WIN32)
  
ENDIF()
SET(FORCE_64BIT_BUILD ON)

SET(boost_GIT_TAG "origin/master")
SET(boost_GIT_URL "https://github.com/CIBC-Internal/boost.git")
SET(BOOST_DIR ${CMAKE_BINARY_DIR}/externals/boost CACHE INTERNAL "boost project directory")

# TODO: fix install step
#
# If CMake ever allows overriding the checkout command or adding flags,
# git checkout -q will silence message about detached head (harmless).
ExternalProject_Add(Boost_external
  GIT_REPOSITORY ${boost_GIT_URL}
  GIT_TAG ${boost_GIT_TAG}
  BUILD_IN_SOURCE ON
  PATCH_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  PREFIX ${BOOST_DIR}
  CMAKE_CACHE_ARGS
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=${CMAKE_VERBOSE_MAKEFILE}
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DFORCE_64BIT_BUILD:BOOL=${FORCE_64BIT_BUILD}
    -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
)

ExternalProject_Get_Property(Boost_external INSTALL_DIR)
ExternalProject_Get_Property(Boost_external SOURCE_DIR)
SET(SCI_BOOST_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "boost include directory")
SET(SCI_BOOST_LIBRARY_DIR ${SOURCE_DIR}/lib)

SET(BOOST_PREFIX "boost_")
SET(THREAD_POSTFIX "-mt")

SET(SCI_BOOST_LIBRARY_TMP) # create a tmp var to append libs to

FOREACH(lib ${boost_Libraries})
  SET(LIB_NAME "${BOOST_PREFIX}${lib}${THREAD_POSTFIX}")
  LIST(APPEND SCI_BOOST_LIBRARY_TMP ${LIB_NAME})
ENDFOREACH()

SET(SCI_BOOST_LIBRARY ${SCI_BOOST_LIBRARY_TMP} CACHE INTERNAL "-Boost library.") # set the lib variable

# Boost is special case - normally this should be handled in external library repo
#CONFIGURE_FILE(${SUPERBUILD_DIR}/BoostConfig.cmake.in ${INSTALL_DIR}/BoostConfig.cmake @ONLY)
#CONFIGURE_FILE(${SUPERBUILD_DIR}/UseBoost.cmake ${SCI_BOOST_USE_FILE} COPYONLY)

SET(Boost_DIR ${INSTALL_DIR} CACHE PATH "")

include_directories(${SCI_BOOST_INCLUDE})
link_directories(${SCI_BOOST_LIBRARY_DIR})
add_definitions(-DBOOST_ALL_NO_LIB)

if(UNIX)
  add_definitions(-DBOOST_NO_CXX11_ALLOCATOR)
endif()


# TODO: if static runtime link is supported, then ABI tag postfix must include s
# see:
# http://www.boost.org/doc/libs/1_58_0/more/getting_started/windows.html
# http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html
if(WIN32)
  set(DEBUG_POSTFIX "-gyd")
  set(boost_LIB_PREFIX "lib")
else()
  set(DEBUG_POSTFIX "-yd")
  set(boost_LIB_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
endif()

foreach(lib ${SCI_BOOST_LIBRARY})
  add_library(${lib} STATIC IMPORTED GLOBAL)

  if(CMAKE_GENERATOR MATCHES "Makefiles" AND CMAKE_BUILD_TYPE MATCHES "Debug")
    set(FULL_LIB_NAME "${SCI_BOOST_LIBRARY_DIR}/${boost_LIB_PREFIX}${lib}${DEBUG_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  elseif(CMAKE_GENERATOR MATCHES "Makefiles")
    set(FULL_LIB_NAME "${SCI_BOOST_LIBRARY_DIR}/${boost_LIB_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  else()
    set(FULL_LIB_NAME_RELEASE "${SCI_BOOST_LIBRARY_DIR}/${boost_LIB_PREFIX}${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    set(FULL_LIB_NAME_DEBUG "${SCI_BOOST_LIBRARY_DIR}/${boost_LIB_PREFIX}${lib}${DEBUG_POSTFIX}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  endif()

  if(CMAKE_GENERATOR MATCHES "Makefiles")
    set_target_properties(${lib}
      PROPERTIES
        IMPORTED_LOCATION ${FULL_LIB_NAME}
    )
  else() # IDEs: Xcode, VS, others...
    set_target_properties(${lib}
      PROPERTIES
        IMPORTED_CONFIGURATIONS "RELEASE;DEBUG"
        IMPORTED_LOCATION_RELEASE ${FULL_LIB_NAME_RELEASE}
        IMPORTED_LOCATION_DEBUG ${FULL_LIB_NAME_DEBUG}
    )
  endif()
endforeach()
