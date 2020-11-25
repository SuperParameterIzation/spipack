#include "spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <MUQ/Optimization/NLoptOptimizer.h>

#include "spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
density(std::make_shared<DensityEstimation>(this->samples, options["DensityOptions"].as<YAML::Node>(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
exponentPara(options["BandwidthExponent"].as<double>(defaults.exponentPara)),
neig(options["NumEigenvalues"].as<std::size_t>(defaults.neig)),
eigensolverTol(options["EigensolverTolerance"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIterations"].as<std::size_t>(defaults.eigensolverMaxIt))
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(std::make_shared<DensityEstimation>(this->samples, options["DensityOptions"].as<YAML::Node>(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
exponentPara(options["BandwidthExponent"].as<double>(defaults.exponentPara)),
neig(options["NumEigenvalues"].as<std::size_t>(defaults.neig)),
eigensolverTol(options["EigensolverTolerance"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIterations"].as<std::size_t>(defaults.eigensolverMaxIt))
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(std::make_shared<DensityEstimation>(this->samples, options["DensityOptions"].as<YAML::Node>(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
exponentPara(options["BandwidthExponent"].as<double>(defaults.exponentPara)),
neig(options["NumEigenvalues"].as<std::size_t>(defaults.neig)),
eigensolverTol(options["EigensolverTolerance"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIterations"].as<std::size_t>(defaults.eigensolverMaxIt))
{}

Eigen::VectorXd KolmogorovOperator::EstimateDensity(bool const tune) const {
  assert(density);

  // estimate the density at each sample
  const Eigen::VectorXd dens = density->Estimate(tune);

  // the tuning (may have) changed the manifold dimension estimation
  manifoldDim = density->ManifoldDimension();

  return dens;
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* tune) {
  assert(kmat.rows()==NumSamples());
  assert(kmat.cols()==NumSamples());

  // estimate the density at each sample
  const Eigen::VectorXd dens = EstimateDensity(*(const bool*)tune);

  // if we are truncating but for some reason want a dense matrix
  if( truncateKernelMatrix ) {
    Eigen::SparseMatrix<double> kmatSparse;
    const Eigen::MatrixXd rowsum = KernelMatrix(eps, dens, kmatSparse);
    kmat = Eigen::MatrixXd(kmatSparse);
    return rowsum;
  }

  return KernelMatrix(eps, dens, kmat);
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::Ref<Eigen::MatrixXd> kmat) {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(dens.rows()==n);

  // compute (and store) the diagonal matrix P
  P = dens.array().pow(exponentPara);

  // compute the unnormalized kernel matrix K
  Eigen::MatrixXd rowsum = SampleRepresentation::KernelMatrix(4.0*eps, P, kmat);

  const double para = VariableBandwidthExponent();

  // compute the normalization for the new kernel matrix
  rowsum = rowsum.array()/P.array().pow(manifoldDim*para);

  // renormalize the kernel matrix
  #pragma omp parallel for num_threads(numThreads)
  for( std::size_t i=0; i<n; ++i ) {
    for( std::size_t j=0; j<n; ++j ) { kmat(i, j) /= rowsum(i)*rowsum(j); }
  }

  return kmat.rowwise().sum();
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* tune) {
  // estimate the density at each sample
  const Eigen::VectorXd dens = EstimateDensity(*(const bool*)tune);

  if( truncateKernelMatrix ) { return KernelMatrix(eps, dens, kmat); }

  // compute a dense kernel matrix using the bandwidth
  Eigen::MatrixXd kmatDense(NumSamples(), NumSamples());
  const Eigen::MatrixXd rowsum = KernelMatrix(eps, dens, kmatDense);

  // convert to sparse form
  kmat = kmatDense.sparseView();
  return rowsum;
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::SparseMatrix<double>& kmat) {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(dens.size()==n);

  const double para = VariableBandwidthExponent();

  // compute (and store) the diagonal matrix P
  P = dens.array().pow(exponentPara);

  // compute the unnormalized kernel matrix K
  std::vector<Eigen::Triplet<double> > entries;
  Eigen::MatrixXd rowsum = SampleRepresentation::KernelMatrix(4.0*eps, P, entries);

  // compute the normalization for the new kernel matrix
  rowsum = rowsum.array()/P.array().pow(manifoldDim*para);

  // renormalize the kernel matrix
  #pragma omp parallel for num_threads(numThreads)
  for( auto& entry : entries ) {
    entry = Eigen::Triplet(entry.row(), entry.col(), entry.value()/(rowsum(entry.row())*rowsum(entry.col())));
  }

  // resize the matrix
  kmat.resize(n, n);
  kmat.setFromTriplets(entries.begin(), entries.end());

  // recompute the rowsum
  rowsum = Eigen::VectorXd::Zero(n);
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }
  return rowsum;
}

void KolmogorovOperator::TuneBandwidthParameter(bool const tuneDens) {
  const Eigen::VectorXd rho = EstimateDensity(tuneDens).array().pow(exponentPara);

  // the cost function and optimizater
  auto cost = std::make_shared<BandwidthCost>(rho, Density());
  auto opt = std::make_shared<NLoptOptimizer>(cost, pt);

  // the initial condition for the optimization is the current parameter value
  std::vector<Eigen::VectorXd> inputs(1);
  assert(bandwidthPara>0.0);
  inputs[0] = Eigen::VectorXd::Constant(1, std::log2(bandwidthPara));

  // solve the optimization and update the parameters
  std::pair<Eigen::VectorXd, double> soln = opt->Solve(inputs);
  bandwidthPara = std::pow(2, soln.first(0))/4.0;
}

double KolmogorovOperator::VariableBandwidthExponent() const { return 1.0+0.5*manifoldDim*exponentPara+exponentPara-0.5*operatorConstant; }

double KolmogorovOperator::ExponentParameter() const { return exponentPara; }

std::shared_ptr<DensityEstimation> KolmogorovOperator::Density() const { return density; }

void KolmogorovOperator::ComputeEigendecomposition(Eigen::Ref<Eigen::VectorXd> Sinv, Eigen::Ref<Eigen::VectorXd> eigenvalues, Eigen::Ref<Eigen::MatrixXd> eigenvectors) {
  // the number of samples
  const std::size_t n = NumSamples();

  assert(eigenvalues.size()==neig);
  assert(eigenvectors.rows()==n);
  assert(eigenvectors.cols()==neig);

  // this function does not tune the density estimation
  const bool tuneDens = false;

  // compute the kernel matrix and the diagonal matrix S^{-1}
  Eigen::SparseMatrix<double> kmat;
  Sinv = KernelMatrix(BandwidthParameter(), kmat, &tuneDens);
  Sinv = (P.array()*Sinv.array().sqrt()).inverse();

  // compute Lhat, which is related to the Kolmogorov operator by a similarity transformation
  Eigen::SparseMatrix<double> Lhat(n, n);
  Lhat = (P.array()*P.array()).inverse().matrix().asDiagonal();
  Lhat = (Sinv.asDiagonal()*kmat*Sinv.asDiagonal()-Lhat)/BandwidthParameter();

  // compute the eigendecomposition of Lhat
  Spectra::SparseGenMatProd<double> matvec(Lhat);
  Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolver(&matvec, neig, std::min(10*neig, n));
  eigsolver.init();
  const std::size_t ncomputed = eigsolver.compute(eigensolverMaxIt, eigensolverTol);
  assert(ncomputed==neig);

  // fill the eigenvectors/values
  eigenvalues = eigsolver.eigenvalues();
  eigenvectors = eigsolver.eigenvectors();
}

std::size_t KolmogorovOperator::NumEigenvalues() const { return neig; }

double KolmogorovOperator::EigensolverTolerance() const { return eigensolverTol; }

double KolmogorovOperator::EigensolverMaxIterations() const { return eigensolverMaxIt; }
