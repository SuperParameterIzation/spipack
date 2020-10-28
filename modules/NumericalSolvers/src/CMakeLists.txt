set(SOURCE
  GraphLaplacian.cpp
)

message(STATUS ${SPIPACK_EXTERNAL_LIBRARIES})

add_library(spiNS SHARED ${SOURCE})
target_link_libraries(spiNS ${SPIPACK_EXTERNAL_LIBRARIES})
