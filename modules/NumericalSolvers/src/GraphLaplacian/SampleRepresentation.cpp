#include "spipack/NumericalSolvers/GraphLaplacian/SampleRepresentation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

SampleRepresentation::SampleRepresentation(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
samples(rv, options["NearestNeighbors"])
{}

SampleRepresentation::SampleRepresentation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
samples(samples, options["NearestNeighbors"])
{}

std::size_t SampleRepresentation::NumSamples() const { return samples.NumSamples(); }

Eigen::Ref<Eigen::VectorXd const> SampleRepresentation::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return samples.Point(i);
}
