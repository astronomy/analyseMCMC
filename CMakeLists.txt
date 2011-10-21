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
  set( CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo Profile."
    FORCE )
endif( NOT CMAKE_BUILD_TYPE )
set( CMAKE_CONFIGURATION_TYPES "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "internal" )


# Project name and language:
project( AnalyseMCMC Fortran )


# Use CPack for packaging:
#include( CPack )
#set( CPACK_INSTALL_PREFIX_PATH ${CMAKE_INSTALL_PREFIX_PATH} )

# Search in the CMake/ directory for CMake modules:
list( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake )


# Various compile/optimisation options that we may want to enable:
include( SetCompileOptions )


# Set directories:
# Place the products in their directories:
get_filename_component( Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME )
set( MODULE_DIRECTORY "${CMAKE_SOURCE_DIR}/mod" )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}" )



# Find dependencies:
find_package( LibSUFR REQUIRED )
find_package( PGPLOT )

set( INCLUDE_FLAGS "-I${LibSUFR_INCLUDES}" )  # will be transferred to CompilerFlags



# Set source files:
include( FileList )

# Set FORTRAN compiler flags:
include( CompilerFlags )




# Code version generator:
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
set_target_properties( analyseMCMC PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )


# Installation targets:
#install( TARGETS "AnalyseMCMC" RUNTIME DESTINATION "bin" )

# Place the executable in the project's root directory:
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )


add_test( AnalyseMCMCruns AnalyseMCMC )

