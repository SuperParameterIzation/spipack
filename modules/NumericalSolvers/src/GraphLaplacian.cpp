#include "spipack/NumericalSolvers/GraphLaplacian.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

GraphLaplacian::GraphLaplacian(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
  cloud(SampleRandomVariable(rv, options["NumSamples"].as<std::size_t>())),
  maxLeaf(options["MaxLeaf"].as<std::size_t>(defaults.maxLeaf)),
  kdtree(cloud.StateDim(), cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeaf))
{}

GraphLaplacian::GraphLaplacian(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options) :
  cloud(samples),
  maxLeaf(options["MaxLeaf"].as<std::size_t>(defaults.maxLeaf)),
  kdtree(cloud.StateDim(), cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeaf))
{}

std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> GraphLaplacian::SampleRandomVariable(std::shared_ptr<RandomVariable> const& rv, std::size_t const n) {
  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  return samples;
}

std::size_t GraphLaplacian::NumSamples() const { return cloud.kdtree_get_point_count(); }

std::size_t GraphLaplacian::KDTreeMaxLeaf() const { return maxLeaf; }

const Eigen::VectorXd& GraphLaplacian::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return cloud.Point(i);
}

void GraphLaplacian::BuildKDTree() {
  // (re-)build the kd-tree
  kdtree.buildIndex();
}

void GraphLaplacian::FindNeighbors(Eigen::VectorXd const& x, double const r, std::vector<std::pair<std::size_t, double> >& neighbors) const {
  assert(x.size()==cloud.StateDim());

  // need to use the squared radius
  const double r2 = r*r;

  // unsorted radius set
  nanoflann::RadiusResultSet<double, std::size_t> resultSet(r2, neighbors); // use squared radius because kd tree metric is the squared euclidean distance

  // find the nearest neighbors---neighbors in a specified radius
  kdtree.findNeighbors(resultSet, x.data(), nanoflann::SearchParams());
}

GraphLaplacian::PointCloud::PointCloud(std::shared_ptr<SampleCollection> const& samples) : samples(samples) {}

std::size_t GraphLaplacian::PointCloud::kdtree_get_point_count() const {
  assert(samples);
  return samples->size();
}

double GraphLaplacian::PointCloud::kdtree_get_pt(std::size_t const p, std::size_t const i) const {
  assert(samples);
  assert(p<samples->size());
  assert(i<samples->at(p)->state[0].size());
  return samples->at(p)->state[0][i];
}

const Eigen::VectorXd& GraphLaplacian::PointCloud::Point(std::size_t const i) const {
  assert(samples);
  assert(i<samples->size());
  return samples->at(i)->state[0];
}

std::size_t GraphLaplacian::PointCloud::StateDim() const {
  assert(samples);
  assert(samples->size()>0);
  return samples->at(0)->state[0].size();
}
