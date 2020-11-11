#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

DensityEstimation::DensityEstimation(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim))
{}

DensityEstimation::DensityEstimation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim))
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
  KernelMatrix(bandwidthPara, squaredBandwidth.array().sqrt(), kmat);

  // compute the volume vector
  const Eigen::VectorXd vol = ((NumSamples()*std::pow(M_PI*bandwidthPara, manifoldDim/2.0))*squaredBandwidth.array().pow(manifoldDim/2.0)).array().inverse();

  // apply the matrix vector product to compute the density estimation
  return kmat*vol;
}
