# =================================
# RUFUS
#
# DataMiningProject/CMakeLists.tx
# =================================


#Set up min CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#Set project name
PROJECT(DMP)

SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

IF("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  MESSAGE(SEND_ERROR "In-source builds are not allowed, please build this in the bin directory!")
ENDIF("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

INCLUDE_DIRECTORIES(core)

# add subfolders
ADD_SUBDIRECTORY(cmakefiles)
ADD_SUBDIRECTORY(core)
ADD_SUBDIRECTORY(tools)
