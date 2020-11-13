#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
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

DensityEstimation::DensityEstimation(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim)),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{}

double DensityEstimation::BandwidthParameter() const { return bandwidthPara; }

Eigen::VectorXd DensityEstimation::Estimate() const {
  // compute the squared bandwidth
  Eigen::VectorXd squaredBandwidth = samples->SquaredBandwidth(numNearestNeighbors);

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
  Eigen::VectorXd squaredBandwidth = samples->SquaredBandwidth(numNearestNeighbors);

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

  // maximum sigma prime so far (and corresponding index)
  double maxsigmaprime = -std::numeric_limits<double>::infinity();
  std::size_t lbest;

  // the best matrix so far
  std::shared_ptr<Eigen::SparseMatrix<double> > kmatbest;

  // loop through the candidate bandwidth parameters
  const Eigen::VectorXd bandwidth = squaredBandwidth.array().sqrt();
  tune.candidateBandwidthParameters.resize(ncandidates);
  tune.kernelAvg.resize(ncandidates);
  tune.logKernelAvgDerivative.resize(ncandidates-1);
  for( std::size_t l=0; l<ncandidates; ++l ) {
    tune.candidateBandwidthParameters(l) = std::exp(tune.bandwidthExponent(l));
    assert(tune.candidateBandwidthParameters(l)>0.0);

    // compute the kernel matrix using this condidate
    auto kmat = std::make_shared<Eigen::SparseMatrix<double> >();
    const Eigen::VectorXd rowsum = KernelMatrix(tune.candidateBandwidthParameters(l), bandwidth, *kmat);
    tune.kernelAvg(l) = rowsum.sum()/n2;

    // determine if this is the best candidate
    const std::size_t lm1 = l-1;
    if( l>0 ) {
      tune.logKernelAvgDerivative(lm1) = (std::log(tune.kernelAvg(l))-std::log(tune.kernelAvg(lm1)))/(std::log(tune.candidateBandwidthParameters(l))-std::log(tune.candidateBandwidthParameters(lm1)));

      if( tune.logKernelAvgDerivative(lm1)>maxsigmaprime ) {
        maxsigmaprime = tune.logKernelAvgDerivative(lm1);
        lbest = lm1;

        kmatbest = kmat;
      }
    }
  }

  // reset the parameters
  bandwidthPara = tune.candidateBandwidthParameters(lbest);
  if( tuneManifoldDimension ) { manifoldDim = 2.0*maxsigmaprime; }

  // compute the volume vector
  const Eigen::VectorXd vol = ((NumSamples()*std::pow(M_PI*bandwidthPara, manifoldDim/2.0))*squaredBandwidth.array().pow(manifoldDim/2.0)).array().inverse();

  // apply the matrix vector product to compute the density estimation
  assert(kmatbest);
  return (*kmatbest)*vol;
}
