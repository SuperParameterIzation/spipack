include(ExternalProject)

set(SPECTRA_DEPENDS )
if( SPIPACK_BUILT_EIGEN3 )
  list(APPEND SPECTRA_DEPENDS EIGEN3)
endif()

ExternalProject_Add(
  SPECTRA
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/spectra
    GIT_REPOSITORY https://github.com/yixuan/spectra.git
    DEPENDS ${SPECTRA_DEPENDS}
    LOG_DOWNLOAD OFF
    LOG_UPDATE OFF
    LOG_PATCH OFF
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
    LOG_TEST OFF
    BUILD_COMMAND ""
  	CONFIGURE_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/external/spectra/include && cp -r ${CMAKE_CURRENT_BINARY_DIR}/external/spectra/src/SPECTRA/include/Spectra ${CMAKE_CURRENT_BINARY_DIR}/external/spectra/include
  )

set(SPECTRA_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/spectra/include/)
