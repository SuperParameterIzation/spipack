#include <gtest/gtest.h>

#include "spipack/Tools/Kernels/ExponentialKernel.hpp"

using namespace spi::Tools;

TEST(ExponentialKernelTests, EvaluateDefault) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;

  // create a hat kernel
  const ExponentialKernel kernel(options);

  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
  x2 = 0.3*x2/x2.norm();
  EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), std::exp(-x2.dot(x2)));
  EXPECT_DOUBLE_EQ(kernel(0.1), std::exp(-0.1));
}

TEST(ExponentialKernelTests, Evaluate) {
  const unsigned int dim = 5;

  // parameters
  const double mag = 2.0;
  const double scale = 3.0;
  const double expon = 1.5;

  // the options for this kernel
  YAML::Node options;
  options["Magnitude"] = mag;
  options["Scale"] = scale;
  options["Exponent"] = expon;

  // create a hat kernel
  const ExponentialKernel kernel(options);

  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
  x2 = 0.3*x2/x2.norm();
  EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), mag*std::exp(-scale*std::pow(x2.dot(x2), expon)));
  EXPECT_DOUBLE_EQ(kernel(0.1), mag*std::exp(-scale*std::pow(0.1, expon)));
}
