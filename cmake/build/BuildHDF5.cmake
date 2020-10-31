include(ExternalProject)

if( NOT DEFINED SPIPACK_INTERNAL_HDF5_VERSION )
  set(SPIPACK_INTERNAL_HDF5_VERSION "1.8.19")
endif()

if(NOT HDF5_EXTERNAL_SOURCE)
  set(HDF5_EXTERNAL_SOURCE https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-${SPIPACK_INTERNAL_HDF5_VERSION}/src/CMake-hdf5-${SPIPACK_INTERNAL_HDF5_VERSION}.tar.gz)
  message(STATUS "Will download HDF5 from ${HDF5_EXTERNAL_SOURCE} during compile.")
endif()

if( SPIPACK_USE_OPENMPI )
	set(HDF5_PARALLEL_FLAG	"--enable-parallel" "CFLAGS=-fPIC -I${MPI_INCLUDE_DIR}")
endif()

ExternalProject_Add(
  HDF5
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/hdf5
    URL ${HDF5_EXTERNAL_SOURCE}
    CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}; ${CMAKE_CURRENT_BINARY_DIR}/external/hdf5/src/HDF5/hdf5-${SPIPACK_INTERNAL_HDF5_VERSION}/configure  CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${HDF5_PARALLEL_FLAG} --prefix=${CMAKE_BINARY_DIR}/external/hdf5 --enable-production --disable-examples
    BUILD_COMMAND make -j5 install &> hdf5-build.txt
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)

set(HDF5_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/external/hdf5/include)
set(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})

set(HDF5_LIBRARY ${CMAKE_BINARY_DIR}/external/hdf5/lib/${library_prefix}hdf5${shared_library_suffix})
set(HDF5HL_LIBRARY ${CMAKE_BINARY_DIR}/external/hdf5/lib/${library_prefix}hdf5_hl${shared_library_suffix})
set(HDF5_LIBRARIES ${HDF5_LIBRARY} ${HDF5HL_LIBRARY})
