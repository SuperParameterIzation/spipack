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
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{
  Initialize(options);
}

DensityEstimation::DensityEstimation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{
  Initialize(options);
}

DensityEstimation::DensityEstimation(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
tuneManifoldDimension(options["TuneManifoldDimension"].as<bool>(defaults.tuneManifoldDimension))
{
  Initialize(options);
}

void DensityEstimation::Initialize(YAML::Node const& options) {
  const YAML::Node& opt = options["Optimization"].as<YAML::Node>(YAML::Node());

  pt.put("Ftol.AbsoluteTolerance", opt["Ftol.AbsoluteTolerance"].as<double>(1.0e-6));
  pt.put("Ftol.RelativeTolerance", opt["Ftol.RelativeTolerance"].as<double>(1.0e-6));
  pt.put("Xtol.AbsoluteTolerance", opt["Xtol.AbsoluteTolerance"].as<double>(1.0e-6));
  pt.put("Xtol.RelativeTolerance", opt["Xtol.RelativeTolerance"].as<double>(1.0e-6));
  pt.put("MaxEvaluations", opt["MaxEvaluations"].as<std::size_t>(1000));
  pt.put("Algorithm", opt["Algorithm"].as<std::string>("COBYLA"));
}

double DensityEstimation::BandwidthParameter() const { return bandwidthPara; }

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

    // the initial condition for the optimization is the current parameter value
    std::vector<Eigen::VectorXd> inputs(1);
    inputs[0] = Eigen::VectorXd::Constant(1, std::log2(bandwidthPara));

    // solve the optimization and update the parameters
    std::pair<Eigen::VectorXd, double> soln = opt->Solve(inputs);
    bandwidthPara = std::pow(2, soln.first(0));
    if( tuneManifoldDimension ) { manifoldDim = -2.0*soln.second; }
  }

  // compute the kernel matrix
  Eigen::SparseMatrix<double> kmat;
  const Eigen::VectorXd rowsum = KernelMatrix(bandwidthPara, bandwidth, kmat);

  // compute the volume vector
  const Eigen::VectorXd vol = ((NumSamples()*std::pow(M_PI*bandwidthPara, manifoldDim/2.0))*squaredBandwidth.array().pow(manifoldDim/2.0)).array().inverse();

  // apply the matrix vector product to compute the density estimation
  return kmat*vol;
}
