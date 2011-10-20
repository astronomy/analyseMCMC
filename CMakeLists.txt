##
##  CMakeFile for AnalyseMCMC
##  MvdS, 24/07/2010
##
##  To compile, from the directory that contains this file, do:
##    $ mkdir build; cd build
##    $ cmake ..
##    $ make
##


cmake_minimum_required( VERSION 2.6 FATAL_ERROR )

# Set build type. Do this *before* we set the project name:
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo Profile."
    FORCE )
endif( NOT CMAKE_BUILD_TYPE )
set( CMAKE_CONFIGURATION_TYPES "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "internal" )


# Project name and language:
project( AnalyseMCMC Fortran )


# Search in the CMake/ directory for CMake modules:
list( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake )


# Various compile/optimisation options that we may want to enable:
include( SetCompileOptions )


# Set directories:
#set( SRCDIR ${CMAKE_SOURCE_DIR}/src )
set( MODDIR ${CMAKE_SOURCE_DIR}/mod )
#set( OBJDIR ${CMAKE_SOURCE_DIR}/obj )
set( BINDIR ${CMAKE_SOURCE_DIR} )


# Define source files:
include( FileList )

# Set search locations:
include( CMakeLocations )


# Find dependencies:
find_package( LibSUFR REQUIRED )
find_package( PGPLOT )

set( INCLUDE_FLAGS "-I${LibSUFR_INCLUDES}" )  # will be transferred to CompilerFlags


# Set FORTRAN compiler flags:
include( CompilerFlags )



# SVN version generator:
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_SOURCE_DIR}/src/code_version.f90
  COMMAND cd $(CMAKE_SOURCE_DIR)
  COMMAND . ${CMAKE_SOURCE_DIR}/code_version.sh $(CMAKE_SOURCE_DIR)/src/code_version.f90 ${Fortran_COMPILER_NAME} ${OPT_FLAGS}
  )

# Tell CMake the source won't be available until build time.
SET_SOURCE_FILES_PROPERTIES( ${CMAKE_SOURCE_DIR}/src/code_version.f90 PROPERTIES GENERATED 1 )



# AnalyseMCMC executable:
add_executable( analyseMCMC ${AnalyseMCMC_SRC_FILES} src/code_version.f90 )
target_link_libraries( analyseMCMC ${LibSUFR_LIBRARIES} ${PGPLOT_LIBRARIES} )
set_target_properties( analyseMCMC PROPERTIES Fortran_MODULE_DIRECTORY ${MODDIR} )


# Installation targets:
#install( TARGETS "AnalyseMCMC" RUNTIME DESTINATION "bin" )

# Place the executable in the project's root directory:
set( EXECUTABLE_OUTPUT_PATH ${BINDIR} )


add_test( AnalyseMCMCruns AnalyseMCMC )

