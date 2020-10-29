include(ExternalProject)

set(MUQ_DEPENDS )

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

list(APPEND MUQ_INCLUDE_DIRS
  "${CMAKE_BINARY_DIR}/external/muq/include"
)

list(APPEND MUQ_LIBRARIES
  "${CMAKE_BINARY_DIR}/external/muq/lib/${library_prefix}muqUtilities${shared_library_suffix}"
  "${CMAKE_BINARY_DIR}/external/muq/lib/${library_prefix}muqModeling${shared_library_suffix}"
  "${CMAKE_BINARY_DIR}/external/muq/lib/${library_prefix}muqSamplingAlgorithms${shared_library_suffix}"
)
