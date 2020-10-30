include(ExternalProject)

set(MUQ_DEPENDS )
if( SPIPACK_BUILT_EIGEN3 )
  list(APPEND SPIPACK_DEPENDS EIGEN3)
endif()
if( SPIPACK_BUILT_HDF5 )
  list(APPEND MUQ_DEPENDS HDF5)
endif()
if( SPIPACK_BUILT_BOOST )
  list(APPEND MUQ_DEPENDS BOOST)
endif()

ExternalProject_Add(
  MUQ
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/muq
    GIT_REPOSITORY https://bitbucket.org/mituq/muq2/src/master/
    DEPENDS ${MUQ_DEPENDS}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/muq
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DMUQ_EIGEN3_DIR=${EIGEN3_INCLUDE_DIR}
    -DMUQ_BOOST_DIR=${SPIPACK_BOOST_DIR}
    -DBOOST_INCLUDE_DIR=${BOOST_INCLUDE_DIRS}
    -DBOOST_SYSTEM_LIBRARY=${BOOST_SYSTEM_LIBRARY}
    -DBOOST_FILESYSTEM_LIBRARY=${BOOST_FILESYSTEM_LIBRARY}
    -DBOOST_GRAPH_LIBRARY=${BOOST_GRAPH_LIBRARY}
    -DMUQ_HDF5_DIR=${SPIPACK_HDF5_DIR}
    -DHDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIR}
    -DHDF5_LIBRARY=${HDF5_LIBRARY}
    -DHDF5HL_LIBRARY=${HDF5HL_LIBRARY}
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
