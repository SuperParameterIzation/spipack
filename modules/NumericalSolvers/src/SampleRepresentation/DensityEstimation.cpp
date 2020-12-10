#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

#include <MUQ/Optimization/NLoptOptimizer.h>

#include "spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

DensityEstimation::DensityEstimation(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{}

DensityEstimation::DensityEstimation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{}

DensityEstimation::DensityEstimation(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{}

Eigen::VectorXd DensityEstimation::Estimate(bool const tune) {
  assert(samples);

  // compute the squared bandwidth
  const Eigen::VectorXd squaredBandwidth = samples->SquaredBandwidth(numNearestNeighbors);

  // compute the density
  return Estimate(squaredBandwidth, tune);
}

Eigen::VectorXd DensityEstimation::Estimate(Eigen::Ref<const Eigen::VectorXd> const& squaredBandwidth, bool const tune) {
  const Eigen::VectorXd bandwidth = squaredBandwidth.array().sqrt();

  if( tune ) {
    // the cost function and optimizater
    auto cost = std::make_shared<BandwidthCost>(bandwidth, shared_from_this());
    auto opt = std::make_shared<NLoptOptimizer>(cost, pt);

    std::cout << "(density) method: " << pt.get("Algorithm", "") << std::endl;

    // the initial condition for the optimization is the current parameter value
    std::vector<Eigen::VectorXd> inputs(1);
    inputs[0] = Eigen::VectorXd::Constant(1, std::log2(bandwidthPara));

    std::cout << "(density) before: " << inputs[0] << std::endl;

    // solve the optimization and update the parameters
    std::pair<Eigen::VectorXd, double> soln = opt->Solve(inputs);
    std::cout << "(density) after: " << soln.first(0) << std::endl;

    bandwidthPara = std::pow(2, soln.first(0));

    if( tuneManifoldDimension ) { manifoldDim = -2.0*soln.second; }
  }

  // compute the kernel matrix
  Eigen::SparseMatrix<double> kmat;
  KernelMatrix(bandwidthPara, bandwidth, kmat);
  assert(kmat.rows()==NumSamples());
  assert(kmat.cols()==NumSamples());

  // compute the volume vector
  Eigen::VectorXd vol = ((NumSamples()*std::pow(M_PI*bandwidthPara, manifoldDim/2.0))*squaredBandwidth.array().pow(manifoldDim/2.0)).array().inverse();
  assert(vol.size()==NumSamples());

  // apply the matrix vector product to compute the density estimation
  return kmat*vol;
}
