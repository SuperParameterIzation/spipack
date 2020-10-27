find_package(MUQ HINTS ${SPIPACK_MUQ_DIR})

if( MUQ_FOUND )
  list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS
    ${MUQ_INCLUDE_DIRS}
  )

  list(APPEND SPIPACK_EXTERNAL_LIBRARIES
    ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES}
  )
else()
  message(STATUS "NO MUQ FOUND")

  include(ExternalProject)

  ExternalProject_Add(MUQ
    GIT_REPOSITORY https://bitbucket.org/mituq/muq2/src/master/
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external
  )
endif()
