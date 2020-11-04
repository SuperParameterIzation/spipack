set(SOURCE
  GraphLaplacian/GraphLaplacian.cpp
)

message(STATUS ${SPIPACK_EXTERNAL_LIBRARIES})

add_library(spiNS SHARED ${SOURCE})
target_link_libraries(spiNS ${SPIPACK_EXTERNAL_LIBRARIES})
add_dependencies(spiNS spiTS)
add_dependencies(spiNS ${SPIPACK_BUILT_DEPENDENCIES})

install(TARGETS spiNS
        EXPORT ${CMAKE_PROJECT_NAME}Depends
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
