##############################################################################
#
# Filename: CMakeLists.txt.h
# Initial Author: Jan Ehrhardt
# Copyright (C) 2010 Institute of Medical Informatics, University of Luebeck
#
##############################################################################

PROJECT(SliceInterpolation)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0012 NEW)
ENDIF(COMMAND CMAKE_POLICY)

##########################################################   
# Set the place where the lib and binary files should be  built

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

##########################################################   
# Set INCLUDE PATHES for THIRD PARTY libraries
#

# TCLAP -- command line parsing
SET(TCLAP_PATH ${PROJECT_SOURCE_DIR}/tclap)
FIND_PATH(TCLAP_INCLUDE_PATH NAMES tclap/CmdLine.h PATHS ${TCLAP_PATH}/include)

##########################################################   
#  Find where ITK was built or installed.

FIND_PACKAGE(ITK)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

SET(ITK_LIBS 
  ${ITK_LIBRARIES}
)        

##########################################################   
#  Find where ITK was built or installed.

INCLUDE_DIRECTORIES(
  ${TCLAP_INCLUDE_PATH}
)

LINK_DIRECTORIES(
  ${LIBRARY_OUTPUT_PATH}
)

# Link ITK libraries
LINK_LIBRARIES(
  ${ITK_LIBS}
)

ADD_EXECUTABLE(imiSliceInterpolation imiSliceInterpolationMain.cxx)


