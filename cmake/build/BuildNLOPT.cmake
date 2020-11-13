include(ExternalProject)

ExternalProject_Add(
  NLOPT
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/nlopt/
    GIT_REPOSITORY https://github.com/stevengj/nlopt.git
    LOG_DOWNLOAD OFF
    LOG_UPDATE OFF
    LOG_PATCH OFF
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
    LOG_TEST OFF
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/external/nlopt
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
)

set(NLOPT_INCLUDE_DIR
  "${CMAKE_BINARY_DIR}/external/nlopt/include"
)

set(NLOPT_LIBRARY
  "${CMAKE_BINARY_DIR}/external/nlopt/lib/${library_prefix}nlopt${shared_library_suffix}"
)
