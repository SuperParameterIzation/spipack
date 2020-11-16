#include "spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
density(std::make_shared<DensityEstimation>(this->samples, options["DensityOptions"].as<YAML::Node>(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
variableBandwidthExponent(options["BandwidthExponent"].as<double>(defaults.variableBandwidthExponent))
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(std::make_shared<DensityEstimation>(this->samples, options["DensityOptions"].as<YAML::Node>(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
variableBandwidthExponent(options["BandwidthExponent"].as<double>(defaults.variableBandwidthExponent))
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(std::make_shared<DensityEstimation>(this->samples, options["DensityOptions"].as<YAML::Node>(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
variableBandwidthExponent(options["BandwidthExponent"].as<double>(defaults.variableBandwidthExponent))
{}

Eigen::VectorXd KolmogorovOperator::EstimateDensity(bool const tune) const {
  assert(density);

  // estimate the density at each sample
  const Eigen::VectorXd dens = density->Estimate(tune);

  // the tuning (may have) changed the manifold dimension estimation
  manifoldDim = density->ManifoldDimension();

  return dens;
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* tune) const {
  // estimate the density at each sample
  const Eigen::VectorXd dens = EstimateDensity(*(const bool*)tune);

  // if we are truncating but for some reason want a dense matrix
  if( truncateKernelMatrix ) {
    Eigen::SparseMatrix<double> kmatSparse;
    const Eigen::VectorXd rowsum = KernelMatrix(eps, dens, kmatSparse);
    kmat = Eigen::MatrixXd(kmatSparse);
    return rowsum;
  }

  return KernelMatrix(eps, dens, kmat);
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::Ref<Eigen::MatrixXd> kmat) const {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(dens.size()==n);

  // compute the unnormalized kernel matrix K
  Eigen::VectorXd rowsum = SampleRepresentation::KernelMatrix(eps, dens.array().pow(variableBandwidthExponent), kmat);

  const double para = 1.0+0.5*manifoldDim*variableBandwidthExponent+variableBandwidthExponent-0.5*operatorConstant;

  // compute the normalization for the new kernel matrix
  rowsum = (rowsum.array()/dens.array().pow(variableBandwidthExponent*manifoldDim)).pow(para);

  // renormalize the kernel matrix
  #pragma omp parallel for num_threads(numThreads)
  for( std::size_t i=0; i<n; ++i ) {
    for( std::size_t j=0; j<n; ++j ) { kmat(i, j) /= rowsum(i)*rowsum(j); }
  }

  return kmat.rowwise().sum();
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* tune) const {
  // estimate the density at each sample
  const Eigen::VectorXd dens = EstimateDensity(*(const bool*)tune);

  if( truncateKernelMatrix ) { return KernelMatrix(eps, dens, kmat); }

  // compute a dense kernel matrix using the bandwidth
  Eigen::MatrixXd kmatDense(NumSamples(), NumSamples());
  const Eigen::VectorXd rowsum = KernelMatrix(eps, dens, kmatDense);

  // convert to sparse form
  kmat = kmatDense.sparseView();
  return rowsum;
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::SparseMatrix<double>& kmat) const {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(dens.size()==n);

  const double para = 1.0+0.5*manifoldDim*variableBandwidthExponent+variableBandwidthExponent-0.5*operatorConstant;

  // compute the unnormalized kernel matrix K
  std::vector<Eigen::Triplet<double> > entries;
  Eigen::VectorXd rowsum = SampleRepresentation::KernelMatrix(eps, dens.array().pow(variableBandwidthExponent), entries);

  // compute the normalization for the new kernel matrix
  rowsum = (rowsum.array()/dens.array().pow(manifoldDim*variableBandwidthExponent)).pow(para);

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
