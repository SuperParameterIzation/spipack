set(SOURCE
  MacroscaleInformation.cpp
  ConditionalVelocityDistribution.cpp

  KineticModels.cpp
  KineticParticleModels.cpp
)

add_library(spiKE SHARED ${SOURCE})
target_link_libraries(spiKE spiNS spiTS ${SPIPACK_EXTERNAL_LIBRARIES})
if( DEFINED SPIPACK_BUILT_DEPENDENCIES )
    add_dependencies(spiKE ${SPIPACK_BUILT_DEPENDENCIES})
endif()

install(TARGETS spiKE
                EXPORT ${CMAKE_PROJECT_NAME}Depends
                LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
                ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
