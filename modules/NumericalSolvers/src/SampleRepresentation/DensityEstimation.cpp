#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

#include <MUQ/Optimization/NLoptOptimizer.h>

#include "spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::Optimization;
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

Eigen::VectorXd DensityEstimation::Estimate(bool const tune) {
  // compute the squared bandwidth
  Eigen::VectorXd squaredBandwidth = samples->SquaredBandwidth(numNearestNeighbors);

  // compute the density
  return Estimate(squaredBandwidth, tune);
}

Eigen::VectorXd DensityEstimation::Estimate(Eigen::Ref<const Eigen::VectorXd> const& squaredBandwidth, bool const tune) {
  const Eigen::VectorXd bandwidth = squaredBandwidth.array().sqrt();

  if( tune ) {
    std::cout << "TUNING THE PARAMETER!" << std::endl;
    // the cost function
    auto cost = std::make_shared<BandwidthCost>(bandwidth, shared_from_this());

    pt::ptree pt;
    pt.put("Optimization.Ftol.AbsoluteTolerance", 1.0e-6);
    pt.put("Optimization.Ftol.RelativeTolerance", 1.0e-6);
    pt.put("Optimization.Xtol.AbsoluteTolerance", 1.0e-6);
    pt.put("Optimization.Xtol.RelativeTolerance", 1.0e-6);
    pt.put("Optimization.MaxEvaluations", 1000); // max number of cost function evaluations
    pt.put("Optimization.Algorithm", "COBYLA");

    auto opt = std::make_shared<NLoptOptimizer>(cost, pt.get_child("Optimization"));

    std::vector<Eigen::VectorXd> inputs(1);
    inputs[0] = Eigen::VectorXd::Constant(1, std::log2(bandwidthPara));

    std::pair<Eigen::VectorXd, double> soln = opt->Solve(inputs);

    std::cout << "opt eps: " << std::pow(2, soln.first(0)) << std::endl;
    std::cout << "cost: " << soln.second << std::endl;

    bandwidthPara = std::pow(2, soln.first(0));
    if( tuneManifoldDimension ) { manifoldDim = 2.0*soln.second; }
  }

  // compute the kernel matrix
  Eigen::SparseMatrix<double> kmat;
  const Eigen::VectorXd rowsum = KernelMatrix(bandwidthPara, squaredBandwidth.array().sqrt(), kmat);

  // compute the volume vector
  const Eigen::VectorXd vol = ((NumSamples()*std::pow(M_PI*bandwidthPara, manifoldDim/2.0))*squaredBandwidth.array().pow(manifoldDim/2.0)).array().inverse();

  // apply the matrix vector product to compute the density estimation
  return kmat*vol;
}
