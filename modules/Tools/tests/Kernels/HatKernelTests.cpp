#include <gtest/gtest.h>

#include "spipack/Tools/Kernels/HatKernel.hpp"

using namespace spi::Tools;

TEST(HatKernelTests, KernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "HatKernel";

  // create a hat kernel
  std::shared_ptr<Kernel> kernel = Kernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check inside the kernel support
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
  x2 = 0.3*x2/x2.norm();
  EXPECT_DOUBLE_EQ(kernel->Evaluate(x1, x2), 1.0);
  EXPECT_DOUBLE_EQ(kernel->operator()(x1, x2), 1.0);
}

TEST(HatKernelTests, IsotropicKernelConstructor) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "HatKernel";

  // create a hat kernel
  std::shared_ptr<IsotropicKernel> kernel = IsotropicKernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check inside the kernel support
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
  x2 = 0.3*x2/x2.norm();
  EXPECT_DOUBLE_EQ(kernel->Evaluate(x1, x2), 1.0);
  EXPECT_DOUBLE_EQ(kernel->operator()(0.1), 1.0);
}

TEST(HatKernelTests, CompactKernelConstructor) {
  const unsigned int dim = 5;

  // the magnitude of the kernel
  const double mag = 2.0;

  // the options for this kernel
  YAML::Node options;
  options["Kernel"] = "HatKernel";
  options["Magnitude"] = mag;

  // create a hat kernel
  std::shared_ptr<CompactKernel> kernel = CompactKernel::Construct(options);

  EXPECT_TRUE(kernel);

  // check inside the kernel support
  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
  x2 = 0.3*x2/x2.norm();  EXPECT_DOUBLE_EQ(kernel->Evaluate(x1, x2), mag);
  EXPECT_DOUBLE_EQ(kernel->operator()(0.1), mag);
}

TEST(HatKernelTests, EvaluateDefault) {
  const unsigned int dim = 5;

  // the options for this kernel
  YAML::Node options;

  // create a hat kernel
  const HatKernel kernel(options);

  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  { // check inside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 0.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), 1.0);
    EXPECT_DOUBLE_EQ(kernel.EvaluateIsotropicKernel(0.1), 1.0);
    EXPECT_DOUBLE_EQ(kernel.EvaluateCompactKernel(0.1), 1.0);
    EXPECT_DOUBLE_EQ(kernel(0.1), 1.0);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1), kernel.IsotropicKernelDerivative(0.1), 1.0e-5);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelSecondDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1), kernel.IsotropicKernelSecondDerivative(0.1), 1.0e-5);
  }

  { // check outside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 1.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), 0.0);
    EXPECT_DOUBLE_EQ(kernel.EvaluateIsotropicKernel(1.1), 0.0);
    EXPECT_DOUBLE_EQ(kernel.EvaluateCompactKernel(1.1), 0.0);
    EXPECT_DOUBLE_EQ(kernel(1.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1), kernel.IsotropicKernelDerivative(0.1), 1.0e-5);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelSecondDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1), kernel.IsotropicKernelSecondDerivative(0.1), 1.0e-5);
  }
}

TEST(HatKernelTests, Evaluate) {
  const unsigned int dim = 5;

  // the magnitude of the kernel
  const double mag = 2.0;

  // the options for this kernel
  YAML::Node options;
  options["Magnitude"] = mag;

  // create a hat kernel
  const HatKernel kernel(options);

  const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
  { // check inside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 0.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), mag);
    EXPECT_DOUBLE_EQ(kernel.EvaluateCompactKernel(0.1), mag);
    EXPECT_DOUBLE_EQ(kernel.EvaluateIsotropicKernel(0.1), mag);
    EXPECT_DOUBLE_EQ(kernel(0.1), mag);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1), kernel.IsotropicKernelDerivative(0.1), 1.0e-5);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelSecondDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1), kernel.IsotropicKernelSecondDerivative(0.1), 1.0e-5);
  }

  { // check outside the kernel support
    Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
    x2 = 1.3*x2/x2.norm();
    EXPECT_DOUBLE_EQ(kernel.Evaluate(x1, x2), 0.0);
    EXPECT_DOUBLE_EQ(kernel.EvaluateIsotropicKernel(1.1), 0.0);
    EXPECT_DOUBLE_EQ(kernel.EvaluateCompactKernel(1.1), 0.0);
    EXPECT_DOUBLE_EQ(kernel(1.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelDerivativeFD(0.1), kernel.IsotropicKernelDerivative(0.1), 1.0e-5);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1),  0.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(kernel.IsotropicKernelSecondDerivative(0.1), 0.0);
    EXPECT_NEAR(kernel.IsotropicKernelSecondDerivativeFD(0.1), kernel.IsotropicKernelSecondDerivative(0.1), 1.0e-5);
  }
}

TEST(HatKernelTests, Serialize) {
  const unsigned int dim = 5;

  // parameters
  const double mag = 2.0;

  std::stringstream buffer;
  {
    // the options for this kernel
    YAML::Node options;
    options["Magnitude"] = mag;

    // create an exponential kernel
    auto kernel = std::make_shared<HatKernel>(options);

    // load kernel to the buffer
    cereal::BinaryOutputArchive oarchive(buffer);
    oarchive(kernel);
  }

  {
    // create an input archive
    cereal::BinaryInputArchive iarchive(buffer);

    // load the kernel from the buffer
    std::shared_ptr<HatKernel> kern;
    iarchive(kern);

    // check kernel parameters
    EXPECT_DOUBLE_EQ(kern->Magnitude(), mag);

    const Eigen::VectorXd x1 = Eigen::VectorXd::Zero(dim, 1);
    { // check inside the kernel support
      Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
      x2 = 0.3*x2/x2.norm();
      EXPECT_DOUBLE_EQ(kern->Evaluate(x1, x2), mag);
      EXPECT_DOUBLE_EQ(kern->EvaluateCompactKernel(0.1), mag);
      EXPECT_DOUBLE_EQ(kern->EvaluateIsotropicKernel(0.1), mag);
      EXPECT_DOUBLE_EQ(kern->operator()(0.1), mag);
    }

    { // check outside the kernel support
      Eigen::VectorXd x2 = Eigen::VectorXd::Random(dim, 1);
      x2 = 1.3*x2/x2.norm();
      EXPECT_DOUBLE_EQ(kern->Evaluate(x1, x2), 0.0);
      EXPECT_DOUBLE_EQ(kern->EvaluateIsotropicKernel(1.1), 0.0);
      EXPECT_DOUBLE_EQ(kern->EvaluateCompactKernel(1.1), 0.0);
      EXPECT_DOUBLE_EQ(kern->operator()(1.1), 0.0);
    }
  }
}
