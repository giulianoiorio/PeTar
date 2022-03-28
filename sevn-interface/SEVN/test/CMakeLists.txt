project(SEVN_test)

include_directories(include include/general include/binstar include/utility)
add_definitions(-DCATCH_CONFIG_ENABLE_BENCHMARKING) #To enable banchmarking

add_executable(test.x sevn_test.cpp)

if(h5)
    target_link_libraries(test.x catch sevn_lib ${HDF5_LIBRARIES})
else()
    target_link_libraries(test.x catch sevn_lib)
endif(h5)


add_custom_command(TARGET test.x POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
        ${CMAKE_SOURCE_DIR}/build/exe/)

add_custom_command(TARGET test.x POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:test.x>
        ${CMAKE_SOURCE_DIR}/build/exe/)
