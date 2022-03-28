
file(GLOB SEVN_BASE "star/*.cpp" "star/lambdas/*.cpp" "star/procs/*.cpp" "general/*.cpp" "general/utils/*.cpp" "binstar/*.cpp" "binstar/procs/*.cpp")
#GI: I removed binstar/procs because it seems that idea is to have all the Processes in Process.h


if(h5)
    file(GLOB H5OUT_BASE "h5out/*.cpp")
    list(APPEND SEVN_BASE ${H5OUT_BASE})
endif(h5)

#message(${H5OUT_BASE})
message(${SEVN_BASE})

add_library(sevn_lib STATIC ${SEVN_BASE})

##install command are enabled with make install
#install(TARGETS sevn_lib DESTINATION /usr/lib)
