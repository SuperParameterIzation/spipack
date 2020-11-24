set(SOURCE
  SampleRepresentation/SampleRepresentation.cpp
  SampleRepresentation/DensityEstimation.cpp
  SampleRepresentation/KolmogorovOperator.cpp

  SampleRepresentation/BandwidthCost.cpp
)

message(STATUS ${SPIPACK_EXTERNAL_LIBRARIES})

add_library(spiNS SHARED ${SOURCE})
target_link_libraries(spiNS spiTS ${SPIPACK_EXTERNAL_LIBRARIES})
if( DEFINED SPIPACK_BUILT_DEPENDENCIES )
    add_dependencies(spiNS spiTS ${SPIPACK_BUILT_DEPENDENCIES})
endif()

install(TARGETS spiNS
        EXPORT ${CMAKE_PROJECT_NAME}Depends
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
