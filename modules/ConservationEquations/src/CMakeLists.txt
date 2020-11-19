set(SOURCE
  #CompressibleEuler.cpp
  Boltzmann.cpp
)

add_library(spiCE SHARED ${SOURCE})
#target_link_libraries(spiCE ${SPIEX_LIBRARIES} ${SPIPACK_EXTERNAL_LIBRARIES})

target_link_libraries(spiCE spiEX_Boltzmann ${SPIPACK_EXTERNAL_LIBRARIES})

if( DEFINED SPIPACK_BUILT_DEPENDENCIES )
  add_dependencies(spiCE spiEX_Boltzmann ${SPIPACK_BUILT_DEPENDENCIES})
endif()

install(TARGETS spiCE
        EXPORT ${CMAKE_PROJECT_NAME}Depends
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
