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

  # Only enable the parts of MUQ that we really want (i.e., MCMC)
  set(MUQ_ENABLEGROUP_DEFAULT OFF CACHE BOOL "MUQ Default Compilegroup status")
  set(MUQ_ENABLEGROUP_UTILITIES_HDF5 ON CACHE BOOL "Enable MUQ HDF5 interface")

  ExternalProject_Add(MUQ
    GIT_REPOSITORY https://bitbucket.org/mituq/muq2/src/master/
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external -DMUQ_ENABLEGROUP_DEFAULT=OFF -DMUQ_ENABLEGROUP_UTILITIES_HDF5=ON
  )
endif()
