include(ExternalProject)

ExternalProject_Add(
  EIGEN3
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/eigen3
    -DCMAKE_CXX_FLAGS=-w --std=c++11
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DBUILD_SHARED_LIBS=ON
  	BUILD_COMMAND make -j5
    INSTALL_COMMAND make -j5 install
)

if(APPLE)
  set(suffix ".dylib")
else()
  set(suffix ".so")
endif()

set(EIGEN3_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/eigen3/include/eigen3/)
