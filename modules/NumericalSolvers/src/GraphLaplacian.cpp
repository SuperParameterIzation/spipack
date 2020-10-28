#include "spipack/NumericalSolvers/GraphLaplacian.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

GraphLaplacian::GraphLaplacian(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) {
  // add random samples into a sample collection
  samples = std::make_shared<SampleCollection>();
  for( size_t i=0; i<options["NumSamples"].as<unsigned int>(); ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }
}

GraphLaplacian::GraphLaplacian(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options) :
  samples(samples)
{}

unsigned int GraphLaplacian::NumSamples() const {
  assert(samples);
  return samples->size();
}
