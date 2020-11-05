include(ExternalProject)

ExternalProject_Add(
  EIGEN3
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3
	GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  LOG_DOWNLOAD OFF
  LOG_UPDATE OFF
  LOG_PATCH OFF
  LOG_CONFIGURE OFF
  LOG_BUILD OFF
  LOG_INSTALL OFF
  LOG_TEST OFF
	BUILD_COMMAND ""
	CONFIGURE_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/Eigen ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/src/EIGEN3/unsupported ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include
)

set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3/include/)
