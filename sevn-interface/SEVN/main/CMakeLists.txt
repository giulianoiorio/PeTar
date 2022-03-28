project(SEVN_app)


if (sin)
    #MAIN program SEVN single evolution
    add_executable(sevn.x sevn.cpp)

    if(h5)
        target_link_libraries(sevn.x sevn_lib ${HDF5_LIBRARIES})
    else()
        target_link_libraries(sevn.x sevn_lib)
    endif(h5)

    add_custom_command(TARGET sevn.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_SOURCE_DIR}/build/exe/)

    add_custom_command(TARGET sevn.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:sevn.x>
            ${CMAKE_SOURCE_DIR}/build/exe/)
endif(sin)



#MAIN program SEVN binary evolution
if (bin)
    add_executable(sevnB.x sevnB.cpp)

    if(h5)
        target_link_libraries(sevnB.x sevn_lib ${HDF5_LIBRARIES})
    else()
        target_link_libraries(sevnB.x sevn_lib)
    endif(h5)

    add_custom_command(TARGET sevnB.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_SOURCE_DIR}/build/exe/)

    add_custom_command(TARGET sevnB.x POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:sevnB.x>
            ${CMAKE_SOURCE_DIR}/build/exe/)
endif(bin)


