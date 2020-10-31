#include <iostream>

// include the google testing header
#include <gtest/gtest.h>

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  const int res = RUN_ALL_TESTS();

  return res;
}
