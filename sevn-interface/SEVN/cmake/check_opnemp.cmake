#Check if OpenMP is present
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
elseif(NOT OPENMP_FOUND)
    message("-- Error: OpenMP not found")
    return()
endif()
