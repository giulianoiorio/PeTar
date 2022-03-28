cmake_minimum_required(VERSION 2.8)
project(SEVN)


#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)


#check flags
include(cmake/check_flags.cmake)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_FLAGS} -W -Wall  -O3")



#OPTIONS AND SETS
option(h5 "Enables HDF5 output files" OFF)
option(doc "Build documentation" OFF)
option(debug "Enable debug output" OFF)
option(test "Enable units test" OFF)
option(bin "Compile sevnB" ON)
option(sin "Compile sevn" ON)
option(verbose "Enable verbose compilation" OFF)

if(verbose)
    set(CMAKE_VERBOSE_MAKEFILE ON)
    message("-- Verbose compilation enabled")
endif(verbose)


#INCLUDE CMAKE FILES
include(cmake/check_opnemp.cmake)
include(cmake/check_h5.cmake)
include(cmake/check_debug.cmake)
include(cmake/generate_doc.cmake)
include(cmake/check_stdfunctions.cmake)

#include header files
include_directories(include include/star include/star/procs include/star/lambdas include/binstar include/binstar/procs include/general include/general/utils include/general/catch2)
include_directories(include/h5out)
#GI: I removed include/binstar/procs because it seems that idea is to have all the Processes in Process.h

#Subdirectory
add_subdirectory(src)
#add_subdirectory(main)


if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
    message("Please consider to switch to CMake 3.8.0")
endif()



if(test AND ${CMAKE_VERSION} VERSION_LESS "3.2.0")
    message(WARNING "You enable the test compilation, but this is possible only with a Cmake version >=3.2, current version is ${CMAKE_VERSION}. The unit test will be disabled.")
    message("-- UNITS TEST compilation disabled")
elseif(test)
    message("-- Units test compilation enabled")
    add_subdirectory(extern/catch)
    add_subdirectory(test)
else()
    message("-- UNITS TEST compilation disabled")
endif()
