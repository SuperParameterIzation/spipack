#include <gtest/gtest.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

class KolmogorovOperatorTests : public::testing::Test {
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
    kernelOptions["Kernel"] = "ExponentialKernel";

    // set the options for the density estimation
    YAML::Node densityOptions;
    densityOptions["NearestNeighbors"] = nnOptions;
    densityOptions["NumNearestNeighbors"] = nneighs;
    densityOptions["KernelOptions"] = kernelOptions;
    densityOptions["ManifoldDimension"] = (double)dim;

    // set the options for the Kolmogorov operator
    options["NearestNeighbors"] = nnOptions;
    options["NumNearestNeighbors"] = nneighs;
    options["KernelOptions"] = kernelOptions;
    options["DensityOptions"] = densityOptions;
    options["BandwidthParameter"] = eps;
    options["BandwidthExponent"] =  bandwidthExponent;
    options["OperatorParameter"] = operatorParameter;
    options["ManifoldDimension"] = (double)dim;
  }

  /// Make sure everything is what we expect
  virtual void TearDown() override {
    EXPECT_EQ(kolOperator->NumNearestNeighbors(), nneighs);
    EXPECT_EQ(kolOperator->NumSamples(), n);
  }

  /// Create the sample representation from samples
  /**
    \return The sample collection used to create the Kolmogorov operator
  */
  inline std::shared_ptr<SampleCollection> CreateFromSamples() {
    // add random samples into a sample collection
    auto samples = std::make_shared<SampleCollection>();
    for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

    // create the Kolmogorov operator
    kolOperator = std::make_shared<KolmogorovOperator>(samples, options);

    // return the samples
    return samples;
  }

  /// The dimension of state spaces
  const unsigned int dim = 1;

  /// The number of samples
  std::size_t n = 1000;

  /// The number of nearest neighbors
  const std::size_t nneighs = 15;

  /// The bandwidth parameter \f$\epsilon\f$
  double eps = 0.15;

  /// The operator parameter \f$c\f$
  const double operatorParameter = 0.75;

  /// The bandwidth exponent \f$\beta\f$
  const double bandwidthExponent = -0.25;

  /// The random variable that lets us sample from the underlying distribution
  std::shared_ptr<RandomVariable> rv;

  /// The options for the graph Laplacian
  YAML::Node options;

  /// The Kolmogorov operator---use a pointer here so we can initalize it as null
  std::shared_ptr<KolmogorovOperator> kolOperator;
};

TEST_F(KolmogorovOperatorTests, RandomVariableConstruction) {
  // create the Kolmogorov operator
  kolOperator = std::make_shared<KolmogorovOperator>(rv, options);
}

TEST_F(KolmogorovOperatorTests, SampleCollectionConstruction) {
  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-kolOperator->Point(i)).norm(), 0.0, 1.0e-10);
  }
}

TEST_F(KolmogorovOperatorTests, NearestNeighborsConstruction) {
  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  auto neighbors = std::make_shared<NearestNeighbors>(samples, options["NearestNeighbors"]);

  // create the sample representation
  kolOperator = std::make_shared<KolmogorovOperator>(neighbors, options);

  // check to make sure the samples match
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR((samples->at(i)->state[0]-kolOperator->Point(i)).norm(), 0.0, 1.0e-10);
  }
}

TEST_F(KolmogorovOperatorTests, UntruncatedKernelMatrix_Dense) {
  options["TruncateKernelMatrix"] = false;

  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // do we want to tune the bandwidth parameter for the density estimation?
  const bool tuneDensity = true;

  // the kernel matrix
  Eigen::MatrixXd kmat(n, n);
  kolOperator->KernelMatrix(eps, kmat, &tuneDensity);

  const Eigen::MatrixXd rowsum = kolOperator->KernelMatrix(eps, kmat, &tuneDensity);
  EXPECT_EQ(rowsum.rows(), n);

  // compute the expected kernel matrix
  Eigen::MatrixXd kernmatExpected(n, n);
  {
    // compute the density to the beta power (do not return the density estimation parameter)
    const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

    auto kern = kolOperator->Kernel();
    #pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];

        kernmatExpected(i, j) = kern->EvaluateIsotropicKernel(diff.dot(diff)/(4.0*eps*rho(i)*rho(j)));
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }

    const double para = 1.0+0.5*dim*bandwidthExponent+bandwidthExponent-0.5*operatorParameter;
    Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();
    rowsumExpected = rowsumExpected.array()/(rho.array().pow((double)dim)).pow(para);

    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        kernmatExpected(i, j) /= rowsumExpected(i)*rowsumExpected(j);
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }
  }
  const Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();

  #pragma omp parallel num_threads(omp_get_max_threads())
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR(rowsum(i), rowsumExpected(i), 1.0e-10);

    double sum = 0.0;
    for( std::size_t j=0; j<n; ++j ) {
      sum += kmat(i, j);
      EXPECT_NEAR(kmat(i, j), kernmatExpected(i, j), 1.0e-10);
    }
    EXPECT_NEAR(sum, rowsumExpected(i), 1.0e-10);
  }
}

TEST_F(KolmogorovOperatorTests, TruncatedKernelMatrix_Dense) {
  const double tol = 5.0e-2;
  options["TruncationTolerance"] = -std::log(tol);

  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // do we want to tune the bandwidth parameter for the density estimation?
  const bool tuneDensity = true;

  // the kernel matrix
  Eigen::MatrixXd kmat(n, n);
  const Eigen::MatrixXd rowsum = kolOperator->KernelMatrix(eps, kmat, &tuneDensity);
  EXPECT_EQ(rowsum.rows(), n);

  // compute the expected kernel matrix
  Eigen::MatrixXd kernmatExpected(n, n);
  {
    // compute the density to the beta power (do not return the density estimation parameter)
    const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

    auto kern = kolOperator->Kernel();
    #pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];
        const double eval = kern->EvaluateIsotropicKernel(diff.dot(diff)/(4.0*eps*rho(i)*rho(j)));

        kernmatExpected(i, j) = eval<tol? 0.0 : eval;
        kernmatExpected(j, i) = eval<tol? 0.0 : kernmatExpected(i, j);
      }
    }

    const double para = 1.0+0.5*dim*bandwidthExponent+bandwidthExponent-0.5*operatorParameter;
    Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();
    rowsumExpected = rowsumExpected.array()/(rho.array().pow((double)dim)).pow(para);

    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        kernmatExpected(i, j) /= rowsumExpected(i)*rowsumExpected(j);
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }
  }
  const Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();

  #pragma omp parallel num_threads(omp_get_max_threads())
  for( std::size_t i=0; i<n; ++i ) {
    EXPECT_NEAR(rowsum(i), rowsumExpected(i), 1.0e-10);

    double sum = 0.0;
    for( std::size_t j=0; j<n; ++j ) {
      sum += kmat(i, j);
      EXPECT_NEAR(kmat(i, j), kernmatExpected(i, j), 1.0e-10);
    }
    EXPECT_NEAR(sum, rowsumExpected(i), 1.0e-10);
  }
}

TEST_F(KolmogorovOperatorTests, UntruncatedKernelMatrix_Sparse) {
  options["TruncateKernelMatrix"] = false;

  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // do we want to tune the bandwidth parameter for the density estimation?
  const bool tuneDensity = true;

  // the kernel matrix
  Eigen::SparseMatrix<double> kmat;
  const Eigen::MatrixXd rowsum = kolOperator->KernelMatrix(eps, kmat, &tuneDensity);
  EXPECT_EQ(rowsum.rows(), n);

  // compute the expected kernel matrix
  Eigen::MatrixXd kernmatExpected(n, n);
  {
    // compute the density to the beta power (do not return the density estimation parameter)
    const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

    auto kern = kolOperator->Kernel();
    #pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];

        kernmatExpected(i, j) = kern->EvaluateIsotropicKernel(diff.dot(diff)/(4.0*eps*rho(i)*rho(j)));
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }

    const double para = 1.0+0.5*dim*bandwidthExponent+bandwidthExponent-0.5*operatorParameter;
    Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();
    rowsumExpected = rowsumExpected.array()/(rho.array().pow((double)dim)).pow(para);

    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        kernmatExpected(i, j) /= rowsumExpected(i)*rowsumExpected(j);
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }
  }
  const Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();

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

TEST_F(KolmogorovOperatorTests, TruncatedKernelMatrix_Sparse) {
  const double tol = 5.0e-2;
  options["TruncationTolerance"] = -std::log(tol);

  // create the Kolmogorov operator from samples
  auto samples = CreateFromSamples();

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // do we want to tune the bandwidth parameter for the density estimation?
  const bool tuneDensity = true;

  // the kernel matrix
  Eigen::SparseMatrix<double> kmat;
  const Eigen::MatrixXd rowsum = kolOperator->KernelMatrix(eps, kmat, &tuneDensity);
  EXPECT_EQ(rowsum.rows(), n);

  // compute the expected kernel matrix
  Eigen::MatrixXd kernmatExpected(n, n);
  {
    // compute the density to the beta power (do not return the density estimation parameter)
    const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

    auto kern = kolOperator->Kernel();
    #pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];
        const double eval = kern->EvaluateIsotropicKernel(diff.dot(diff)/(4.0*eps*rho(i)*rho(j)));

        kernmatExpected(i, j) = eval<tol? 0.0 : eval;
        kernmatExpected(j, i) = eval<tol? 0.0 : kernmatExpected(i, j);
      }
    }

    const double para = 1.0+0.5*dim*bandwidthExponent+bandwidthExponent-0.5*operatorParameter;
    Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();
    rowsumExpected = rowsumExpected.array()/(rho.array().pow((double)dim)).pow(para);

    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        kernmatExpected(i, j) /= rowsumExpected(i)*rowsumExpected(j);
        kernmatExpected(j, i) = kernmatExpected(i, j);
      }
    }
  }
  const Eigen::VectorXd rowsumExpected = kernmatExpected.rowwise().sum();

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

TEST_F(KolmogorovOperatorTests, KernelDerivativeAverage_Untruncated) {
  options["TruncateKernelMatrix"] = false;

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // compute the scaled density
  const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

  // compute the kernel matrix
  const double eps = 10.0;
  const double avg = kolOperator->KernelDerivativeAverage(eps, rho);

  // compute the expected kernel derivative average
  double expectedAvg = 0.0;
  {
    auto kern = kolOperator->Kernel();
    //#pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];
        const double theta = diff.dot(diff)/(eps*rho(i)*rho(j));
        const double deriv = -kern->IsotropicKernelDerivative(theta)*theta/eps;

        expectedAvg += (i==j? 1.0 : 2.0)*deriv;
        assert(!std::isnan(expectedAvg));
      }
    }
    expectedAvg /= n*n;
  }

  EXPECT_NEAR(expectedAvg, avg, 1.0e-10);
}

TEST_F(KolmogorovOperatorTests, KernelDerivativeAverage_Truncated) {
  const double tol = 5.0e-2;
  options["TruncationTolerance"] = -std::log(tol);

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // compute the scaled density
  const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

  // compute the kernel matrix
  const double eps = 10.0;
  const double avg = kolOperator->KernelDerivativeAverage(eps, rho);

  // compute the expected kernel derivative average
  double expectedAvg = 0.0;
  {
    auto kern = kolOperator->Kernel();
    //#pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];
        const double theta = diff.dot(diff)/(eps*rho(i)*rho(j));
        const double eval = kern->EvaluateIsotropicKernel(theta);
        const double deriv = -kern->IsotropicKernelDerivative(theta)*theta/eps;

        expectedAvg += (eval<tol? 0.0 : (i==j? 1.0 : 2.0)*deriv);
        assert(!std::isnan(expectedAvg));
      }
    }
    expectedAvg /= n*n;
  }

  EXPECT_NEAR(expectedAvg, avg, 1.0e-10);
}

TEST_F(KolmogorovOperatorTests, KernelSecondDerivativeAverage_Untruncated) {
  options["TruncateKernelMatrix"] = false;

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // compute the scaled density
  const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

  // compute the kernel matrix
  const double eps = 10.0;
  const double avg = kolOperator->KernelSecondDerivativeAverage(eps, rho);

  // compute the expected kernel derivative average
  double expectedAvg = 0.0;
  {
    auto kern = kolOperator->Kernel();
    //#pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];
        const double theta = diff.dot(diff)/(eps*rho(i)*rho(j));
        const double deriv = 2.0*kern->IsotropicKernelDerivative(theta)*theta/eps/eps + kern->IsotropicKernelSecondDerivative(theta)*theta*theta/eps/eps ;

        expectedAvg += (i==j? 1.0 : 2.0)*deriv;
        assert(!std::isnan(expectedAvg));
      }
    }
    expectedAvg /= n*n;
  }

  EXPECT_NEAR(expectedAvg, avg, 1.0e-10);
}

TEST_F(KolmogorovOperatorTests, KernelSecondDerivativeAverage_Truncated) {
  const double tol = 5.0e-2;
  options["TruncationTolerance"] = -std::log(tol);

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // compute the scaled density
  const Eigen::VectorXd rho = kolOperator->EstimateDensity(false).array().pow(bandwidthExponent);

  // compute the kernel matrix
  const double eps = 10.0;
  const double avg = kolOperator->KernelSecondDerivativeAverage(eps, rho);

  // compute the expected kernel derivative average
  double expectedAvg = 0.0;
  {
    auto kern = kolOperator->Kernel();
    //#pragma omp parallel num_threads(omp_get_max_threads())
    for( std::size_t i=0; i<n; ++i ) {
      for( std::size_t j=i; j<n; ++j ) {
        const Eigen::VectorXd diff = samples->at(i)->state[0]-samples->at(j)->state[0];
        const double theta = diff.dot(diff)/(eps*rho(i)*rho(j));
        const double eval = kern->EvaluateIsotropicKernel(theta);
        const double deriv = 2.0*kern->IsotropicKernelDerivative(theta)*theta/eps/eps + kern->IsotropicKernelSecondDerivative(theta)*theta*theta/eps/eps ;

        expectedAvg += (eval<tol? 0.0 : (i==j? 1.0 : 2.0)*deriv);
        assert(!std::isnan(expectedAvg));
      }
    }
    expectedAvg /= n*n;
  }

  EXPECT_NEAR(expectedAvg, avg, 1.0e-10);
}

TEST_F(KolmogorovOperatorTests, Tuning) {
  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  kolOperator->TuneBandwidthParameter();
  EXPECT_TRUE(kolOperator->BandwidthParameter()>0.0);
}

TEST_F(KolmogorovOperatorTests, Eigendecomposition) {
  options["NumEigenvalues"] = 15;
  options["EigensolverTolerance"] = 1.0e-8;
  options["EigensolverMaxIterations"] = 1e5;

  // create the graph laplacian from samples
  auto samples = CreateFromSamples();
  EXPECT_EQ(samples->size(), n);

  // check the eigendecomposition parameters
  EXPECT_EQ(kolOperator->NumEigenvalues(), 15);
  EXPECT_DOUBLE_EQ(kolOperator->EigensolverTolerance(), 1.0e-8);
  EXPECT_EQ(kolOperator->EigensolverMaxIterations(), 1e5);

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // compute the eigendecomposition
  Eigen::VectorXd Sinv(kolOperator->NumSamples());
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::Random(kolOperator->NumEigenvalues()); // initialize to random so we check to make sure the small (magnitude) is set to zero
  Eigen::MatrixXd eigenvectors(kolOperator->NumSamples(), kolOperator->NumEigenvalues());
  kolOperator->ComputeEigendecomposition(Sinv, eigenvalues, eigenvectors);

  // the smallest eigenvalue is zero
  EXPECT_NEAR(eigenvalues(0), 0.0, 1.0e-8);
}
