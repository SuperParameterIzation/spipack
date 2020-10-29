#include "spipack/NumericalSolvers/GraphLaplacian.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

GraphLaplacian::GraphLaplacian(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
  cloud(SampleRandomVariable(rv, options["NumSamples"].as<size_t>())),
  maxLeaf(options["MaxLeaf"].as<size_t>(defaults.maxLeaf))
{}

GraphLaplacian::GraphLaplacian(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options) :
  cloud(samples),
  maxLeaf(options["MaxLeaf"].as<size_t>(defaults.maxLeaf))
{}

std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> GraphLaplacian::SampleRandomVariable(std::shared_ptr<RandomVariable> const& rv, size_t const n) {
  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  return samples;
}

size_t GraphLaplacian::NumSamples() const { return cloud.kdtree_get_point_count(); }

size_t GraphLaplacian::KDTreeMaxLeaf() const { return maxLeaf; }

GraphLaplacian::PointCloud::PointCloud(std::shared_ptr<SampleCollection> const& samples) : samples(samples) {}

size_t GraphLaplacian::PointCloud::kdtree_get_point_count() const {
  assert(samples);
  return samples->size();
}

double GraphLaplacian::PointCloud::kdtree_get_pt(size_t const p, size_t const i) const {
  assert(samples);
  assert(p<samples->size());
  assert(i<samples->at(p)->state[0].size());
  return samples->at(p)->state[0][i];
}
