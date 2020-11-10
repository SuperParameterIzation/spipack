#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/UniformBox.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/Tools/Kernels/HatKernel.hpp"

#include "spipack/NumericalSolvers/SampleRepresentation/GraphLaplacian.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

class GraphLaplacianTests : public::testing::Test {
public:
  /// Set up information to test the graph Laplacian
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // options for the nearest neighbor search
    YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;

    // set the options for the graph laplacian
    options["NearestNeighbors"] = nnOptions;
    options["BandwidthRange.Min"] = bandwidthRange.first;
    options["BandwidthRange.Max"] = bandwidthRange.second;
    options["NumBandwidthSteps"] = numBandwidthSteps;
    options["BandwidthIndex"] = bandwidthIndex;
    options["EigensolverTol"] = eigensolverTol;

    // set the kernel options
    YAML::Node kernelOptions;
    kernelOptions["Kernel"] = "HatKernel";
    options["KernelOptions"] = kernelOptions;
  }

  /// Make sure everything is constructed correctly
  virtual void TearDown() override {
    EXPECT_TRUE(laplacian);

    // make sure the graph laplacian has enough samples
    EXPECT_EQ(laplacian->NumSamples(), n);

    // check the kernel
    EXPECT_TRUE(laplacian->Kernel());
    std::shared_ptr<HatKernel const> kern = std::dynamic_pointer_cast<HatKernel const>(laplacian->Kernel());
    EXPECT_TRUE(kern);

    // check the bandwidth
    EXPECT_NEAR(laplacian->BandwidthRange().first, bandwidthRange.first, 1.0e-10);
    EXPECT_NEAR(laplacian->BandwidthRange().second, bandwidthRange.second, 1.0e-10);
    EXPECT_EQ(laplacian->NumBandwidthSteps(), numBandwidthSteps);


    const Eigen::VectorXd candidateBandwidths = laplacian->BandwidthParameterCandidates();
    EXPECT_EQ(candidateBandwidths.size(), numBandwidthSteps+1);
    EXPECT_NEAR(candidateBandwidths(0), std::pow(2.0, bandwidthRange.first), 1.0e-10);
    EXPECT_NEAR(candidateBandwidths(numBandwidthSteps), std::pow(2.0, bandwidthRange.second), 1.0e-10);
  }

protected:

  /// Create the graph Laplacian from samples
  /**
    \return The sample collection used to create the graph Laplacian
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the graph laplacian
    laplacian = std::make_shared<GraphLaplacian>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const unsigned int dim = 4;

  /// The number of samples
  const std::size_t n = 2000;

  /// The range for the bandwidth parameter \f$2^{l}\f$ (this is the range of \f$l\f$)
  const std::pair<double, double> bandwidthRange = std::pair<double, double>(-20.0, 10.0);

  /// The number of steps in the discretization of the bandwidth parameter range
  const std::size_t numBandwidthSteps = 25;

  // The bandwidth index parameter
  const int bandwidthIndex = 4;

  /// The tolerance for the eigensolver
  const double eigensolverTol = 1.0e-6;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The graph Laplacian---use a pointer here so we can initalize it as null
  std::shared_ptr<GraphLaplacian> laplacian;
};

TEST_F(GraphLaplacianTests, RandomVariableConstruction) {
  // create the graph laplacian
  laplacian = std::make_shared<GraphLaplacian>(rv, options);
}

TEST_F(GraphLaplacianTests, SampleCollectionConstruction) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-laplacian->Point(i)).norm(), 0.0, 1.0e-10);
  }
}

TEST_F(GraphLaplacianTests, ConstructKernelMatrix) {
  // the kernel bandwidth parameter
  const double eps = 0.75;

  // the number of nearest neighbors---used to define the variable bandwidth
  const std::size_t numNeighbors = 10;

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // build the kd tree
  laplacian->BuildKDTrees();

  // compute the bandwidth
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  const Eigen::VectorXd squaredBandwidth = laplacian->SquaredBandwidth(neighbors);

  // compute the kernel matrix
  Eigen::SparseMatrix<double> kernmat(samples->size(), samples->size());
  EXPECT_EQ(kernmat.nonZeros(), 0);
  laplacian->KernelMatrix(eps, squaredBandwidth.array().sqrt(), kernmat);
  EXPECT_TRUE(kernmat.nonZeros()<samples->size()*samples->size());
  EXPECT_TRUE(kernmat.nonZeros()>=samples->size());

  // compute the approximate kernel matrix
  Eigen::SparseMatrix<double> kernmatApprox(samples->size(), samples->size());
  EXPECT_EQ(kernmatApprox.nonZeros(), 0);
  laplacian->KernelMatrix(eps, squaredBandwidth.array().sqrt(), neighbors, kernmatApprox);
  EXPECT_TRUE(kernmatApprox.nonZeros()<samples->size()*samples->size());
  EXPECT_TRUE(kernmatApprox.nonZeros()>=samples->size());

  // compute the expected kernel matrix
  Eigen::MatrixXd kernmatExpected(samples->size(), samples->size());
  {
    auto kern = laplacian->Kernel();
    for( std::size_t i=0; i<samples->size(); ++i ) {
      for( std::size_t j=i; j<samples->size(); ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];

        kernmatExpected(i, j) = kern->EvaluateCompactKernel(diff.dot(diff)/(eps*std::sqrt(squaredBandwidth(i)*squaredBandwidth(j))));
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }
  }

  for( std::size_t i=0; i<samples->size(); ++i ) {
    for( std::size_t j=0; j<samples->size(); ++j ) {
      EXPECT_NEAR(kernmatExpected(i, j), kernmatExpected(j, i), 1.0e-10);
      EXPECT_NEAR(kernmat.coeff(i, j), kernmat.coeff(j, i), 1.0e-10);
      EXPECT_NEAR(kernmatApprox.coeff(i, j), kernmatApprox.coeff(j, i), 1.0e-10);
      EXPECT_NEAR(kernmatExpected(i, j), kernmat.coeff(i, j), 1.0e-10);
      if( kernmatApprox.coeff(i, j)>0.0 ) {
        EXPECT_NEAR(kernmatExpected(i, j), kernmatApprox.coeff(i, j), 1.0e-10);
      } else {
        EXPECT_NEAR(0.0, kernmatApprox.coeff(i, j), 1.0e-10);
      }
    }
  }
}

TEST_F(GraphLaplacianTests, DensityEstimation_TunedBandwidth) {
  options["TuneBandwidth"] = true;

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // check the bandwidth
  EXPECT_TRUE(laplacian->TuneBandwidthParameter());

  // build the kd tree
  laplacian->BuildKDTrees();

  // estimate the density
  const Eigen::VectorXd density = laplacian->DensityEstimation();
}

TEST_F(GraphLaplacianTests, DensityEstimation_FixedBandwidth) {
  options["TuneBandwidth"] = false;

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // check the bandwidth
  EXPECT_FALSE(laplacian->TuneBandwidthParameter());

  // build the kd tree
  laplacian->BuildKDTrees();

  // estimate the density
  const Eigen::VectorXd density = laplacian->DensityEstimation();
}


/*TEST_F(GraphLaplacianTests, ConstructHeatMatrix) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // construct the heat matrix
  laplacian->ConstructHeatMatrix();

  // the sum of the heat matrix rows should be 1
  const Eigen::Ref<const Eigen::SparseMatrix<double> > heat = laplacian->HeatMatrix();
  for( std::size_t i=0; i<laplacian->NumSamples(); ++i ) {
    double sum = 0.0;
    for( std::size_t j=0; j<laplacian->NumSamples(); ++j ) {
      sum += heat.coeff(i, j);
    }
    EXPECT_NEAR(sum, 1.0, 1.0e-10);
  }

  // the largest eigenvalue is 1
  const std::size_t neig = 1;
  const Eigen::VectorXd eigenvalues = laplacian->HeatMatrixEigenvalues(neig);
  EXPECT_NEAR(eigenvalues(0), 1.0, 10.0*eigensolverTol);
}*/

TEST(WeightedPoissonProblem, Solve) {
  /*// the dimension of the problem
  const std::size_t dim = 6;

  // create a random variable
  std::vector<std::pair<double, double> > bounds(dim, std::pair<double, double>(1.0, 0.0));
  auto rv = std::make_shared<UniformBox>(bounds)->AsVariable();

  // the options for the graph Laplacian
  YAML::Node options;
  options["NumSamples"] = 1000;
  options["Bandwidth"] = 0.25;
  //options["EigensolverMaxIt"] = 10000;
  options["EigensolverTol"] = 1.0e-4;

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "HatKernel";
  options["KernelOptions"] = kernelOptions;

  // create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(rv, options);

  // create the right hand side for the weighted poisson problem
  Eigen::VectorXd rhs = Eigen::VectorXd::Ones(laplacian->NumSamples());

  // solve the weighted Poisson problem
  laplacian->SolveWeightedPoisson(rhs);

  std::cout << rhs(0) << std::endl;*/
}
