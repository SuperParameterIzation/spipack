#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

DensityEstimation::DensityEstimation(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim)),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{}

DensityEstimation::DensityEstimation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim)),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{}

double DensityEstimation::BandwidthParameter() const { return bandwidthPara; }

Eigen::VectorXd DensityEstimation::Estimate() const {
  // compute the squared bandwidth
  Eigen::VectorXd squaredBandwidth = samples.SquaredBandwidth(numNearestNeighbors);

  // compute the density
  return Estimate(squaredBandwidth);
}

Eigen::VectorXd DensityEstimation::Estimate(Eigen::Ref<const Eigen::VectorXd> const& squaredBandwidth) const {
  // compute the kernel matrix
  Eigen::SparseMatrix<double> kmat;
  const Eigen::VectorXd rowsum = KernelMatrix(bandwidthPara, squaredBandwidth.array().sqrt(), kmat);

  // compute the volume vector
  const Eigen::VectorXd vol = ((NumSamples()*std::pow(M_PI*bandwidthPara, manifoldDim/2.0))*squaredBandwidth.array().pow(manifoldDim/2.0)).array().inverse();

  // apply the matrix vector product to compute the density estimation
  return kmat*vol;
}

Eigen::VectorXd DensityEstimation::Estimate(DensityEstimation::TuningData& tune) {
  // compute the squared bandwidth
  Eigen::VectorXd squaredBandwidth = samples.SquaredBandwidth(numNearestNeighbors);

  // compute the density
  return Estimate(squaredBandwidth, tune);
}

Eigen::VectorXd DensityEstimation::Estimate(Eigen::Ref<const Eigen::VectorXd> const& squaredBandwidth, DensityEstimation::TuningData& tune) {
  // the number of candidates
  const std::size_t ncandidates = tune.bandwidthExponent.size();

  // the number of samples squared
  const std::size_t n2 = NumSamples()*NumSamples();

  // if there are no candidate bandwidth parameters, just use the fixed parameter
  if( ncandidates==0 ) { return Estimate(squaredBandwidth); }

  // loop through the candidate bandwidth parameters
  const Eigen::VectorXd bandwidth = squaredBandwidth.array().sqrt();
  tune.candidateBandwidthParameters.resize(ncandidates);
  tune.kernelAvg.resize(ncandidates);
  tune.logKernelAvgDerivative.resize(ncandidates-1);
  for( std::size_t l=0; l<ncandidates; ++l ) {
    tune.candidateBandwidthParameters(l) = std::exp(tune.bandwidthExponent(l));
    assert(tune.candidateBandwidthParameters(l)>0.0);

    // compute the kernel matrix using this condidate
    Eigen::SparseMatrix<double> kmat;
    const Eigen::VectorXd rowsum = KernelMatrix(tune.candidateBandwidthParameters(l), bandwidth, kmat);
    tune.kernelAvg(l) = rowsum.sum()/n2;
    if( l>0 ) { tune.logKernelAvgDerivative(l-1) = (std::log(tune.kernelAvg(l))-std::log(tune.kernelAvg(l-1)))/(std::log(tune.candidateBandwidthParameters(l))-std::log(tune.candidateBandwidthParameters(l-1))); }
  }

  // reset the parameters
  std::size_t maxind;
  const double maxsigprime = tune.logKernelAvgDerivative.maxCoeff(&maxind);
  bandwidthPara = tune.candidateBandwidthParameters(maxind);
  if( tuneManifoldDimension ) { manifoldDim = 2.0*maxsigprime; }

  return Estimate(squaredBandwidth);
}
