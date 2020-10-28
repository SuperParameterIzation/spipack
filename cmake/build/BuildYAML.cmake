include(ExternalProject)

ExternalProject_Add(
  YAML
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/yaml
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/yaml
    -DCMAKE_INSTALL_LIBDIR=${CMAKE_BINARY_DIR}/external/yaml/lib
    -DCMAKE_CXX_FLAGS=-w --std=c++11
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DBUILD_SHARED_LIBS=ON
    -DYAML_CPP_BUILD_TESTS=OFF
  	BUILD_COMMAND make
    INSTALL_COMMAND make install
)

if(APPLE)
  set(suffix ".dylib")
else()
  set(suffix ".so")
endif()

set(YamlCpp_LIBRARY ${CMAKE_BINARY_DIR}/external/yaml/lib/libyaml-cpp${suffix})
set(YamlCpp_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/yaml/include)
