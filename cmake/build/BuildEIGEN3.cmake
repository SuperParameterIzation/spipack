include(ExternalProject)

set(EIGEN_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3")

if(NOT EIGEN_EXTERNAL_SOURCE)
	#set(EIGEN_EXTERNAL_SOURCE https://gitlab.com/libeigen/eigen/-/archive/3.3.4.tar.bz2)
	set(EIGEN_EXTERNAL_SOURCE https://gitlab.com/libeigen/eigen/-/archive/3.3.8.tar.bz2)
endif()

set(EIGEN3_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/muq_external/)
ExternalProject_Add(
  EIGEN3
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3
  URL ${EIGEN_EXTERNAL_SOURCE}
	BUILD_COMMAND ""
	CONFIGURE_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/external/eigen3/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/Eigen ${CMAKE_BINARY_DIR}/external/eigen3/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/unsupported ${CMAKE_BINARY_DIR}/external/eigen3/include
)


if(APPLE)
  set(suffix ".dylib")
else()
  set(suffix ".so")
endif()

set(EIGEN3_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/eigen3/include/)
