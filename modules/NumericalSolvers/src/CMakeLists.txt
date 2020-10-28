set(SOURCE
  GraphLaplacian.cpp
)

message(STATUS ${SPIPACK_EXTERNAL_LIBRARIES})

add_library(spiNS SHARED ${SOURCE})
target_link_libraries(spiNS ${SPIPACK_EXTERNAL_LIBRARIES})
if( SPIPACK_BUILT_YAML )
  add_dependencies(spiNS YAML)
endif()
if( SPIPACK_BUILT_MUQ )
  add_dependencies(spiNS MUQ)
endif()
