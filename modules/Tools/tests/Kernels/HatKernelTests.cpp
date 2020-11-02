#include <gtest/gtest.h>

#include "spipack/Tools/Kernels/HatKernel.hpp"

using namespace spi::Tools;

TEST(HatKernelTests, Evaluate) {
  const unsigned int dim = 5;

  // create a hat kernel
  const HatKernel kernel;

  Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);

  { // check inside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
    x2 = 0.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), 1.0);
    EXPECT_DOUBLE_EQ(kernel(0.1), 1.0);
  }

  Eigen::VectorXd diff = x1;

  { // check outside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
    x2 = 1.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), 0.0);
    EXPECT_DOUBLE_EQ(kernel(1.1), 0.0);
  }
}
