set(SOURCE
  GraphLaplacian/GraphLaplacian.cpp
  GraphLaplacian/GraphLaplacianKernel.cpp
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
if( SPIPACK_BUILT_HDF5 )
  add_dependencies(spiNS HDF5)
endif()
if( SPIPACK_BUILT_BOOST )
  add_dependencies(spiNS BOOST)
endif()
if( SPIPACK_BUILT_EIGEN3 )
  add_dependencies(spiNS HDF5)
endif()
