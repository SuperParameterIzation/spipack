include(ExternalProject)

set(MUQ_DEPENDS )
if( SPIPACK_BUILT_EIGEN3 )
  set(MUQ_DEPENDS EIGEN3)
endif()

message(STATUS "")
message(STATUS "")
message(STATUS ${EIGEN3_INCLUDE_DIR})
message(STATUS "")
message(STATUS "")

ExternalProject_Add(
  MUQ
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/muq
    GIT_REPOSITORY https://bitbucket.org/mituq/muq2/src/master/
    DEPENDS ${MUQ_DEPENDS}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/muq
    -DMUQ_EIGEN3_DIR=${EIGEN3_INCLUDE_DIR}
    -DMUQ_ENABLEGROUP_DEFAULT=OFF
    -DMUQ_ENABLEGROUP_UTILITIES_HDF5=ON
    -DMUQ_ENABLEGROUP_SAMPLING_ALGORITHM=ON
  	BUILD_COMMAND make -j5
    INSTALL_COMMAND make -j5 install
)

if(APPLE)
  set(suffix ".dylib")
else()
  set(suffix ".so")
endif()

list(APPEND MUQ_INCLUDE_DIRS
  "${CMAKE_BINARY_DIR}/external/muq/include"
)

list(APPEND MUQ_LIBRARIES
  "${CMAKE_BINARY_DIR}/external/muq/lib/libmuqUtilities${suffix}"
  "${CMAKE_BINARY_DIR}/external/muq/lib/libmuqModeling${suffix}"
  "${CMAKE_BINARY_DIR}/external/muq/lib/libmuqSamplingAlgorithms${suffix}"
)
