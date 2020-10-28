if( NOT DEFINED SPIPACK_GTEST_DIR )
  include(ExternalProject)

  ExternalProject_Add(GTEST
    GIT_REPOSITORY https://github.com/google/googletest.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/gtest
    -DCMAKE_CXX_FLAGS=-Wall --std=c++11
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DBUILD_GMOCK=OFF
  )

  list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS
    "${CMAKE_BINARY_DIR}/external/gtest/include"
  )

  list(APPEND SPIPACK_EXTERNAL_LIBRARIES
    "${CMAKE_BINARY_DIR}/external/gtest/lib64/libgtest.a"
  )
else()
  find_path(GTEST_INCLUDE_DIR gtest/gtest.h
	         HINTS ${SPIPACK_GTEST_DIR}/include ${SPIPACK_GTEST_DIR}/include
	         PATH_SUFFIXES gtest NO_DEFAULT_PATH)

	find_library(GTEST_LIBRARY_STATIC NAMES libgtest.a
	             HINTS ${SPIPACK_GTEST_DIR}/lib ${SPIPACK_GTEST_DIR}/build NO_DEFAULT_PATH)

  list(APPEND SPIPACK_EXTERNAL_INCLUDE_DIRS ${GTEST_INCLUDE_DIR})

  list(APPEND SPIPACK_EXTERNAL_LIBRARIES ${GTEST_LIBRARY_STATIC})
endif()
