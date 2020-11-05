include(ExternalProject)

ExternalProject_Add(
  YAML
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/yaml
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    LOG_DOWNLOAD OFF
    LOG_UPDATE OFF
    LOG_PATCH OFF
    LOG_CONFIGURE OFF
    LOG_BUILD OFF
    LOG_INSTALL OFF
    LOG_TEST OFF
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/yaml
    -DCMAKE_INSTALL_LIBDIR=${CMAKE_BINARY_DIR}/external/yaml/lib
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DBUILD_SHARED_LIBS=ON
    -DYAML_CPP_BUILD_TESTS=OFF
  	BUILD_COMMAND make -j5 install
    INSTALL_COMMAND ""
)

set(YamlCpp_LIBRARY ${CMAKE_BINARY_DIR}/external/yaml/lib/${library_prefix}yaml-cpp${shared_library_suffix})
set(YamlCpp_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/yaml/include)
