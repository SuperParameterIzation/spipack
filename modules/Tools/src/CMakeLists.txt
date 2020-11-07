set(SOURCE
  Kernels/Kernel.cpp
  Kernels/IsotropicKernel.cpp
  Kernels/CompactKernel.cpp
  Kernels/ExponentialKernel.cpp
  Kernels/HatKernel.cpp
  Kernels/BumpKernel.cpp
)

message(STATUS ${SPIPACK_EXTERNAL_LIBRARIES})

add_library(spiTS SHARED ${SOURCE})
target_link_libraries(spiTS ${SPIPACK_EXTERNAL_LIBRARIES})
if( DEFINED SPIPACK_BUILT_DEPENDENCIES )
    add_dependencies(spiTS ${SPIPACK_BUILT_DEPENDENCIES})
endif()

install(TARGETS spiTS
                EXPORT ${CMAKE_PROJECT_NAME}Depends
                LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
                ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
