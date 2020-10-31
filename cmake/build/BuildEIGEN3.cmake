include(ExternalProject)

if(NOT EIGEN_EXTERNAL_SOURCE)
	set(EIGEN_EXTERNAL_SOURCE https://gitlab.com/libeigen/eigen/-/archive/3.3.8.tar.bz2)
endif()

ExternalProject_Add(
  EIGEN3
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3
  URL ${EIGEN_EXTERNAL_SOURCE}
	BUILD_COMMAND ""
	CONFIGURE_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/Eigen ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/unsupported ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include
)

set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include/)
