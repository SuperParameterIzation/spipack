set(SOURCE
  Kernels/Kernel.cpp
  Kernels/IsotropicKernel.cpp
  Kernels/CompactKernel.cpp
  Kernels/HatKernel.cpp
)

message(STATUS ${SPIPACK_EXTERNAL_LIBRARIES})

add_library(spiTS SHARED ${SOURCE})
target_link_libraries(spiTS ${SPIPACK_EXTERNAL_LIBRARIES})
add_dependencies(spiTS ${SPIPACK_BUILT_DEPENDENCIES})

install(TARGETS spiTS
                EXPORT ${CMAKE_PROJECT_NAME}Depends
                LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
                ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
