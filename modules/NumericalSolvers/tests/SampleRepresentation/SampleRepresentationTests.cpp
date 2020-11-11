#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

class SampleRepresentationTests : public::testing::Test {
protected:

  /// Set up information to test the sample representation
  virtual void SetUp() override {
    // create a standard Gaussian random variable
    rv = std::make_shared<Gaussian>(dim)->AsVariable();

    // options for the nearest neighbor search
    YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;
    nnOptions["NumThreads"] = omp_get_max_threads();

    // set the kernel options
    YAML::Node kernelOptions;
    kernelOptions["Kernel"] = "HatKernel";

    // set the options for the graph laplacian
    options["NearestNeighbors"] = nnOptions;
    options["NumNearestNeighbors"] = nneighs;
    options["KernelOptions"] = kernelOptions;
  }

  /// Make sure everything is what we expect
  virtual void TearDown() override {
    EXPECT_EQ(representation->NumNearestNeighbors(), nneighs);
    EXPECT_EQ(representation->NumSamples(), n);
  }

  /// Create the sample representation from samples
  /**
    \return The sample collection used to create the graph Laplacian
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the graph laplacian
    representation = std::make_shared<SampleRepresentation>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const unsigned int dim = 4;

  /// The number of samples
  std::size_t n = 1000;

  /// The number of nearest neighbors
  const std::size_t nneighs = 15;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The sample representation---use a pointer here so we can initalize it as null
  std::shared_ptr<SampleRepresentation> representation;
};

TEST_F(SampleRepresentationTests, RandomVariableConstruction) {
  // create the graph laplacian
  representation = std::make_shared<SampleRepresentation>(rv, options);
}

TEST_F(SampleRepresentationTests, SampleCollectionConstruction) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-representation->Point(i)).norm(), 0.0, 1.0e-10);
  }
}

TEST_F(SampleRepresentationTests, KernelMatrix) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // construct the kd-trees
  representation->BuildKDTrees();

  // compute the kernel matrix
  const double eps = 10.0;
  Eigen::SparseMatrix<double> kmat;
  const Eigen::VectorXd rowsum = representation->KernelMatrix(eps, kmat);
  EXPECT_EQ(rowsum.size(), n);
  EXPECT_EQ(kmat.rows(), n);
  EXPECT_EQ(kmat.cols(), n);
  EXPECT_TRUE(kmat.nonZeros()>=n);
  EXPECT_TRUE(kmat.nonZeros()<n*n);

  // compute the expected kernel matrix
  Eigen::MatrixXd kernmatExpected(n, n);
  {
    // create a nearest neighbors object
    NearestNeighbors nn(samples, options["NearestNeighbors"]);
    std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
    nn.BuildKDTrees();
    const Eigen::VectorXd squaredBandwidth = nn.SquaredBandwidth(nneighs, neighbors);

    auto kern = representation->Kernel();
    #pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];

        kernmatExpected(i, j) = kern->EvaluateIsotropicKernel(diff.dot(diff)/(eps*std::sqrt(squaredBandwidth(i)*squaredBandwidth(j))));
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }
  }
  const Eigen::VectorXd rowsumExpected = kernmatExpected.colwise().sum();

  #pragma omp parallel num_threads(omp_get_max_threads())
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR(rowsum(i), rowsumExpected(i), 1.0e-10);

    double sum = 0.0;
    for( std::size_t j=0; j<n; ++j ) {
      sum += kmat.coeff(i, j);
      EXPECT_NEAR(kmat.coeff(i, j), kernmatExpected(i, j), 1.0e-10);
    }
    EXPECT_NEAR(sum, rowsumExpected(i), 1.0e-10);
  }
}

/*TEST_F(SampleRepresentationTests, DensityEstimation) {
  // reset the number of samples
  n = 10000;

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // compute the kernel matrix
  const double eps = 10.0;
  Eigen::SparseMatrix<double> kmat;
  const Eigen::VectorXd rowsum = representation->KernelMatrix(eps, kmat);
  EXPECT_EQ(rowsum.size(), n);
  EXPECT_EQ(kmat.rows(), n);
  EXPECT_EQ(kmat.cols(), n);
  EXPECT_TRUE(kmat.nonZeros()>=n);
  EXPECT_TRUE(kmat.nonZeros()<n*n);
}*/
