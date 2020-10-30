include(ExternalProject)

ExternalProject_Add(GTEST
  GIT_REPOSITORY https://github.com/google/googletest.git
  STEP_TARGETS install
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/gtest
  -DCMAKE_INSTALL_LIBDIR=${CMAKE_BINARY_DIR}/external/gtest/lib
  -DCMAKE_CXX_FLAGS=-Wall --std=c++11
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  BUILD_COMMAND make -j5 &> gtest-build.txt
  INSTALL_COMMAND make -j5 install
)

set(GTEST_INCLUDE_DIR
  ${CMAKE_BINARY_DIR}/external/gtest/include)

set(GTEST_LIBRARIES
  ${CMAKE_BINARY_DIR}/external/gtest/lib/${library_prefix}gmock${static_library_suffix}
  ${CMAKE_BINARY_DIR}/external/gtest/lib/${library_prefix}gtest${static_library_suffix})
