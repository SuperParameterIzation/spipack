#include "spipack/NumericalSolvers/GraphLaplacian/SampleRepresentation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

SampleRepresentation::SampleRepresentation(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) {}

SampleRepresentation::SampleRepresentation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) {}
