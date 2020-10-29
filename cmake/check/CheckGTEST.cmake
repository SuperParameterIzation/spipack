# make sure that Eigen supports the "Ref" command
set(CMAKE_REQUIRED_INCLUDES ${GTEST_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${GTEST_LIBRARIES})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <gtest/gtest.h>
  TEST(Test, Testing) { EXPECT_TRUE(true); }
  int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
  }
  "
  GTEST_CODE_COMPILES)

set(GTEST_COMPILES 1)
if( NOT GTEST_CODE_COMPILES )
  set(GTEST_COMPILES 0)
endif()
