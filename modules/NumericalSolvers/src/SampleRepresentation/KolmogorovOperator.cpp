#include "spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
density(this->samples, options["DensityOptions"].as<YAML::Node>(options)) // default to using the same parameters as the Kolmogorov operator
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(this->samples, options["DensityOptions"].as<YAML::Node>(options)) // default to using the same parameters as the Kolmogorov operator
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(this->samples, options["DensityOptions"].as<YAML::Node>(options)) // default to using the same parameters as the Kolmogorov operator
{}

void KolmogorovOperator::TuneDensityEstimation() {
  // create the tuning data
  DensityEstimation::TuningData tune;
  tune.bandwidthExponent = Eigen::VectorXd::LinSpaced(15, -5.0, 5.0);

  // estimate the density at each sample
  const Eigen::VectorXd densityEstimate = density.Estimate(tune);

  std::cout << "density estimate: " << densityEstimate.transpose() << std::endl;
}