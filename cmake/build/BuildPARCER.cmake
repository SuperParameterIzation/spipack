include(ExternalProject)

ExternalProject_Add(
  PARCER
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/parcer/
    GIT_REPOSITORY https://andrew_d_davis@bitbucket.org/mituq/parcer.git
    LOG_DOWNLOAD OFF
    LOG_UPDATE OFF
    LOG_PATCH OFF
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
    LOG_TEST OFF
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/external/parcer
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
)

set(PARCER_INCLUDE_DIR
  "${CMAKE_BINARY_DIR}/external/parcer/include"
)

set(PARCER_LIBRARY
  "${CMAKE_BINARY_DIR}/external/parcer/lib/${library_prefix}parcer${shared_library_suffix}"
)
