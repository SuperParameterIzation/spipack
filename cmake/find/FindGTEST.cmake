find_path(GTEST_INCLUDE_DIR gtest/gtest.h
         HINTS ${SPIPACK_GTEST_DIR}/include ${SPIPACK_GTEST_DIR}/include
         PATH_SUFFIXES gtest NO_DEFAULT_PATH)

find_library(GTEST_LIBRARIES NAMES libgtest.a
             HINTS ${SPIPACK_GTEST_DIR}/lib ${SPIPACK_GTEST_DIR}/build NO_DEFAULT_PATH)

if( GTEST_LIBRARIES )
  set(GTEST_FOUND 1)
else()
  set(GTEST_FOUND 0)
endif()
