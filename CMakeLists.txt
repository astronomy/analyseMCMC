##
##  CMakeFile for AnalyseMCMC
##  MvdS, 2010-07-24
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
# libSUFR:
find_package( LibSUFR REQUIRED )
set( INCLUDE_FLAGS "-I${LibSUFR_INCLUDES}" )  # will be transferred to CompilerFlags

# PGPlot:
find_package( PGPLOT )

# Find PLplot if PGPlot isn't found:
if( NOT PGPLOT_FOUND OR FORCE_PLPLOT)
  find_package( PLplot REQUIRED )
  set( INCLUDE_FLAGS "${INCLUDE_FLAGS} -I${PLplot_INCLUDES}" )  # will be transferred to CompilerFlags
endif( NOT PGPLOT_FOUND OR FORCE_PLPLOT )





# Set source files:
include( FileList )

# Set FORTRAN compiler flags:
include( CompilerFlags_Fortran )




# Code-version generator: always create both files:
# This implies that the files code_version* are removed by 'make clean'!
if( CREATE_VERSION )
  
  set_source_files_properties( ${CMAKE_SOURCE_DIR}/src/code_version_pgplot.f90 PROPERTIES GENERATED 1 )  # Tell CMake this source file won't be available until build time.
  set_source_files_properties( ${CMAKE_SOURCE_DIR}/src/code_version_plplot.f90 PROPERTIES GENERATED 1 )  # Tell CMake this source file won't be available until build time.
  
  add_custom_command(
    # PGPlot:
    OUTPUT ${CMAKE_SOURCE_DIR}/src/code_version_pgplot.f90
    COMMAND cd $(CMAKE_SOURCE_DIR)
    COMMAND . ${CMAKE_SOURCE_DIR}/code_version.sh ${CMAKE_SOURCE_DIR} src/code_version_pgplot.f90  no  ${Fortran_COMPILER_NAME} ${OPT_FLAGS}
    
    # PLplot:
    OUTPUT ${CMAKE_SOURCE_DIR}/src/code_version_plplot.f90
    COMMAND cd $(CMAKE_SOURCE_DIR)
    COMMAND . ${CMAKE_SOURCE_DIR}/code_version.sh $(CMAKE_SOURCE_DIR) src/code_version_plplot.f90  yes  ${Fortran_COMPILER_NAME} ${OPT_FLAGS}
    )
  
endif( CREATE_VERSION )


# Use PLplot if found, otherwise use PGPLOT:
if( PLplot_FOUND )
  set( Plot_LIBRARIES ${PLplot_LIBRARIES} )
  message( STATUS "Linking against PLplot" )
else( PLplot_FOUND )
  set( Plot_LIBRARIES ${PGPLOT_LIBRARIES} )
  message( STATUS "Linking against PGPlot" )
endif( PLplot_FOUND )


# Build, link and install AnalyseMCMC executable:
add_executable( analyseMCMC ${AnalyseMCMC_SRC_FILES} ${Plot_SRC_FILES} )
target_link_libraries( analyseMCMC ${LibSUFR_LIBRARIES} ${Plot_LIBRARIES} )
set_target_properties( analyseMCMC PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )
install( TARGETS "analyseMCMC" RUNTIME DESTINATION "bin" )


# Use add_dependencies below in order to pretend other executables depend on aM.  This
#   seems to be necessary to ensure that they are built *after* aM.  If this doesn't
#   happen, several threads will be pouring data into the code_version files
#   simultaneously.

# Build, link and install compPDFs executable:
add_executable( compPDFs ${compPDFs_SRC_FILES} ${Plot_SRC_FILES} )
target_link_libraries( compPDFs ${LibSUFR_LIBRARIES} ${Plot_LIBRARIES} )
set_target_properties( compPDFs PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )
add_dependencies( compPDFs analyseMCMC )  # Ensure it is build after aM
install( TARGETS "compPDFs" RUNTIME DESTINATION "bin" )

# Build, link and install MCMCstats executable:
add_executable( MCMCstats ${MCMCstats_SRC_FILES} ${Plot_SRC_FILES} )
target_link_libraries( MCMCstats ${LibSUFR_LIBRARIES} ${Plot_LIBRARIES} )
set_target_properties( MCMCstats PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )
add_dependencies( MCMCstats analyseMCMC )  # Ensure it is build after aM
install( TARGETS "MCMCstats" RUNTIME DESTINATION "bin" )

# Build, link and install plotSignal executable:
add_executable( plotSignal ${plotSignal_SRC_FILES} ${Plot_SRC_FILES} )
target_link_libraries( plotSignal ${LibSUFR_LIBRARIES} ${Plot_LIBRARIES} )
set_target_properties( plotSignal PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )
add_dependencies( plotSignal analyseMCMC )  # Ensure it is build after aM
install( TARGETS "plotSignal" RUNTIME DESTINATION "bin" )

# Build, link and install plotSpectrum executable:
add_executable( plotSpectrum ${plotSpectrum_SRC_FILES} ${Plot_SRC_FILES} )
target_link_libraries( plotSpectrum ${LibSUFR_LIBRARIES} ${Plot_LIBRARIES} )
set_target_properties( plotSpectrum PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY} )
add_dependencies( plotSpectrum analyseMCMC )  # Ensure it is build after aM
install( TARGETS "plotSpectrum" RUNTIME DESTINATION "bin" )




# Place the executable in the project's root directory:
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )


add_test( AnalyseMCMCruns AnalyseMCMC )

