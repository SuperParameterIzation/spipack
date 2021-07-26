set(SOURCE
  #CompressibleEuler.cpp
  Boltzmann.cpp
)

add_library(spiCE SHARED ${SOURCE})
#target_link_libraries(spiCE ${SPIEX_LIBRARIES} ${SPIPACK_EXTERNAL_LIBRARIES})

target_link_libraries(spiCE spiEX_Boltzmann spiKE spiNS spiTS ${SPIPACK_EXTERNAL_LIBRARIES})
add_dependencies(spiCE spiEX_Boltzmann spiKE spiNS spiTS ${SPIPACK_BUILT_DEPENDENCIES})

if( DEFINED SPIPACK_BUILT_DEPENDENCIES )
  add_dependencies(spiCE ${SPIPACK_BUILT_DEPENDENCIES})
endif()

install(TARGETS spiCE
        EXPORT ${CMAKE_PROJECT_NAME}Depends
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
