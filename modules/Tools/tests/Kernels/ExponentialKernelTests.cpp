#include <gtest/gtest.h>

#include <cereal/archives/binary.hpp>

#include "spipack/Tools/Kernels/ExponentialKernel.hpp"

using namespace spi::Tools;

TEST(ExponentialKernelTests, KernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "ExponentialKernel";

  // create an exponential kernel
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

  // create an exponential kernel
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

  // create an exponential kernel
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

  // create an exponential kernel
  const ExponentialKernel kernel(options);

  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
  x2 = 0.3*x2/x2.norm();
  EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), mag*std::exp(-scale*std::pow(x2.dot(x2), expon)));
  EXPECT_DOUBLE_EQ(kernel.EvaluateIsotropicKernel(0.1), mag*std::exp(-scale*std::pow(0.1, expon)));
  EXPECT_DOUBLE_EQ(kernel(0.1), mag*std::exp(-scale*std::pow(0.1, expon)));
}

TEST(ExponentialKernelTests, Serialize) {
  const unsigned int dim = 5;

  // parameters
  const double mag = 2.0;
  const double scale = 3.0;
  const double expon = 1.5;

  std::stringstream buffer;
  {
    // the options for this kernel
    YAML::Node options;
    options["Magnitude"] = mag;
    options["Scale"] = scale;
    options["Exponent"] = expon;

    // create an exponential kernel
    auto kernel = std::make_shared<ExponentialKernel>(options);

    // load kernel to the buffer
    cereal::BinaryOutputArchive oarchive(buffer);
    oarchive(kernel);
  }

  {
    // create an input archive
    cereal::BinaryInputArchive iarchive(buffer);

    // load the kernel from the buffer
    std::shared_ptr<ExponentialKernel> kern;
    iarchive(kern);

    // check kernel parameters
    EXPECT_DOUBLE_EQ(kern->Magnitude(), mag);
    EXPECT_DOUBLE_EQ(kern->Scale(), scale);
    EXPECT_DOUBLE_EQ(kern->Exponent(), expon);

    // check kernel evaluate
    const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
    Eigen::VectorXd x2 = Eigen::VectorXd::Ones(dim, 1);
    x2 = 0.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kern->Evaluate(x1, x2), mag*std::exp(-scale*std::pow(x2.dot(x2), expon)));
    EXPECT_DOUBLE_EQ(kern->EvaluateIsotropicKernel(0.1), mag*std::exp(-scale*std::pow(0.1, expon)));
    EXPECT_DOUBLE_EQ(kern->operator()(0.1), mag*std::exp(-scale*std::pow(0.1, expon)));
  }
}
