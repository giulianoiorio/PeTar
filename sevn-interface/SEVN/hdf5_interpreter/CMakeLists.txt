cmake_minimum_required(VERSION 3.2)
project(h5inter)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -O3 -std=c++14")


    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
    elseif(NOT HDF5_FOUND)
        message("-- Error: HDF5 not found")
        return()
    endif()


    find_package(Readline)
if(READLINE_FOUND)
	include_directories(${Readline_INCLUDE_DIRS})
        elseif(NOT READLINE_FOUND)
    message("-- Error: Readline library not found")
    return()
endif()

   find_package(Curses REQUIRED)
if(CURSES_FOUND)
    include_directories(${CURSES_INCLUDE_DIR})
elseif(NOT CURSES_FOUND)
    message("-- Error: Curses library not found")
    return()
endif()


find_package(PythonLibs)
if(PYTHONLIBS_FOUND)
    include_directories(${PYTHON_INCLUDE_PATH})
elseif(NOT PYTHONLIBS_FOUND)
    message("-- Error: Python library not found")
    return()
endif()


include_directories(include python_cxx)

file(GLOB H5base "src/*.cpp")
add_library(h5base STATIC ${H5base})


add_executable(run.x ${CMAKE_SOURCE_DIR}/main.cpp)
target_link_libraries(run.x h5base ${HDF5_LIBRARIES} ${Readline_LIBRARIES} ${CURSES_LIBRARIES} ${PYTHON_LIBRARIES})
