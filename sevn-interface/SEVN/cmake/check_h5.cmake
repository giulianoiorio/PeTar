#Enable HDF5 output files if requested
option(h5 "Enables HDF5 output files" ON)
if(h5)
    message("-- HDF5 output files enabled")
    add_definitions(-DH5OUT)

    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
    elseif(NOT HDF5_FOUND)
        message("-- Error: HDF5 not found")
        return()
    endif()

endif(h5)
