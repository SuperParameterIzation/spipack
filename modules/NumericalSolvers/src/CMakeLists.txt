set(SOURCE
  GraphLaplacian/GraphLaplacian.cpp
)

message(STATUS ${SPIPACK_EXTERNAL_LIBRARIES})

add_library(spiNS SHARED ${SOURCE})
target_link_libraries(spiNS ${SPIPACK_EXTERNAL_LIBRARIES})
add_dependencies(spiNS spiTS)
add_dependencies(spiNS ${SPIPACK_BUILT_DEPENDENCIES})
