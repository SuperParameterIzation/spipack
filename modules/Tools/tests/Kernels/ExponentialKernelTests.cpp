#include <gtest/gtest.h>

#include "spipack/Tools/Kernels/ExponentialKernel.hpp"

using namespace spi::Tools;

TEST(ExponentialKernelTests, KernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "ExponentialKernel";

  // create a hat kernel
  std::shared_ptr<Kernel> kernel = Kernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check the kernel
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
  x2 = 0.3*x2/x2.norm();
  EXPECT_DOUBLE_EQ(kernel->Evaluate(x1, x2), std::exp(-x2.dot(x2)));
  EXPECT_DOUBLE_EQ(kernel->operator()(x1, x2), std::exp(-x2.dot(x2)));
}

TEST(ExponentialKernelTests, IsotropicKernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "ExponentialKernel";

  // create a hat kernel
  std::shared_ptr<IsotropicKernel> kernel = IsotropicKernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check the kernel
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
  x2 = 0.3*x2/x2.norm();
  EXPECT_DOUBLE_EQ(kernel->Evaluate(x1, x2), std::exp(-x2.dot(x2)));
  EXPECT_DOUBLE_EQ(kernel->operator()(0.1), std::exp(-0.1));
}

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
  EXPECT_DOUBLE_EQ(kernel.EvaluateIsotropicKernel(0.1), std::exp(-0.1));
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
  EXPECT_DOUBLE_EQ(kernel.EvaluateIsotropicKernel(0.1), mag*std::exp(-scale*std::pow(0.1, expon)));
  EXPECT_DOUBLE_EQ(kernel(0.1), mag*std::exp(-scale*std::pow(0.1, expon)));
}
