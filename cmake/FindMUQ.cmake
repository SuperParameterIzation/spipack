find_package(MUQ HINTS ${SPIPACK_MUQ_DIR})

if( MUQ_FOUND )
  list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS
    ${MUQ_INCLUDE_DIRS}
  )

  list(APPEND SPIPACK_EXTERNAL_LIBRARIES
    ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES}
  )
else()
  include(ExternalProject)

  ExternalProject_Add(MUQ
    GIT_REPOSITORY https://bitbucket.org/mituq/muq2/src/master/
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/muq
    #-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    #-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DMUQ_ENABLEGROUP_DEFAULT=OFF
    -DMUQ_ENABLEGROUP_UTILITIES_HDF5=ON
    -DMUQ_ENABLEGROUP_SAMPLING_ALGORITHM=ON
  )

  if(APPLE)
    set(suffix ".dylib")
  else()
    set(suffix ".so")
  endif()

  list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS
    "${CMAKE_BINARY_DIR}/external/muq/include"
  )

  list(APPEND SPIPACK_EXTERNAL_LIBRARIES
    "${CMAKE_BINARY_DIR}/external/muq/lib/libmuqUtilities${suffix}"
  )
endif()
