#include <gtest/gtest.h>

#include "spipack/Tools/Kernels/BumpKernel.hpp"

using namespace spi::Tools;

TEST(BumpKernelTests, KernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "BumpKernel";

  // create a Bump kernel
  std::shared_ptr<Kernel> kernel = Kernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check inside the kernel support
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
  x2 = 0.3*x2/x2.norm();
  const double theta = x2.dot(x2);
  EXPECT_NEAR(kernel->Evaluate(x1, x2), std::exp(1.0-1.0/(1-theta)), 1.0e-10);
  EXPECT_NEAR(kernel->operator()(x1, x2), std::exp(1.0-1.0/(1-theta)), 1.0e-10);
}

TEST(BumpKernelTests, IsotropicKernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "BumpKernel";

  // create a Bump kernel
  std::shared_ptr<IsotropicKernel> kernel = IsotropicKernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check inside the kernel support
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
  x2 = 0.3*x2/x2.norm();
  const double theta = x2.dot(x2);
  EXPECT_NEAR(kernel->Evaluate(x1, x2), std::exp(1.0-1.0/(1-theta)), 1.0e-10);
  EXPECT_NEAR(kernel->operator()(0.1), std::exp(1.0-1.0/0.9), 1.0e-10);
}

TEST(BumpKernelTests, CompactKernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "BumpKernel";

  // create a bump kernel
  std::shared_ptr<CompactKernel> kernel = CompactKernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check inside the kernel support
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
  x2 = 0.3*x2/x2.norm();
  const double theta = x2.dot(x2);
  EXPECT_NEAR(kernel->Evaluate(x1, x2), std::exp(1.0-1.0/(1-theta)), 1.0e-10);
  EXPECT_NEAR(kernel->operator()(0.1), std::exp(1.0-1.0/0.9), 1.0e-10);
}

TEST(BumpKernelTests, EvaluateDefault) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;

  // create a bump kernel
  const BumpKernel kernel(options);

  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  { // check inside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 0.3*x2/x2.norm();
    const double theta = x2.dot(x2);
    EXPECT_NEAR(kernel.Evaluate(x1, x2), std::exp(1.0-1.0/(1-theta)), 1.0e-10);
    EXPECT_NEAR(kernel(0.1), std::exp(1.0-1.0/0.9), 1.0e-10);
    EXPECT_DOUBLE_EQ(kernel(1.0), 0.0);
  }

  { // check outside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 1.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), 0.0);
    EXPECT_DOUBLE_EQ(kernel(1.1), 0.0);
  }
}

TEST(BumpKernelTests, Evaluate) {
  const unsigned int dim = 5;

  const double mag = 2.0;
  const double scale = 2.5;
  const double expon = 1.25;

  // the options for this kernel
  YAML::Node options;
  options["Magnitude"] = mag;
  options["Scale"] = scale;
  options["Exponent"] = expon;

  // create a bump kernel
  const BumpKernel kernel(options);

  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  { // check inside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 0.3*x2/x2.norm();
    const double theta = x2.dot(x2);
    EXPECT_NEAR(kernel.Evaluate(x1, x2), mag*std::exp(scale*(1.0-1.0/(1-std::pow(theta, expon)))), 1.0e-10);
    EXPECT_NEAR(kernel(0.1), mag*std::exp(scale*(1.0-1.0/(1.0-std::pow(0.1, expon)))), 1.0e-10);
  }

  { // check outside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 1.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), 0.0);
    EXPECT_DOUBLE_EQ(kernel(1.1), 0.0);
  }
}
