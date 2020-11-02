#include "spipack/NumericalSolvers/GraphLaplacian/GraphLaplacian.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

GraphLaplacian::GraphLaplacian(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
  cloud(SampleRandomVariable(rv, options["NumSamples"].as<std::size_t>())),
  maxLeaf(options["MaxLeaf"].as<std::size_t>(defaults.maxLeaf)),
  kdtree(cloud.StateDim(), cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeaf))
{
  Initialize(options);
}

GraphLaplacian::GraphLaplacian(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options) :
  cloud(samples),
  maxLeaf(options["MaxLeaf"].as<std::size_t>(defaults.maxLeaf)),
  kdtree(cloud.StateDim(), cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeaf))
{
  Initialize(options);
}

void GraphLaplacian::Initialize(YAML::Node const& options) {
  // create the kernel
  kernel = CompactKernel::Construct(options["KernelOptions"]);
  assert(kernel);

  // initialize the sparse heat matrix
  heatMatrix.resize(NumSamples(), NumSamples());
}

std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> GraphLaplacian::SampleRandomVariable(std::shared_ptr<RandomVariable> const& rv, std::size_t const n) {
  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  return samples;
}

std::size_t GraphLaplacian::NumSamples() const { return cloud.kdtree_get_point_count(); }

std::size_t GraphLaplacian::KDTreeMaxLeaf() const { return maxLeaf; }

Eigen::Ref<const Eigen::VectorXd> GraphLaplacian::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return cloud.Point(i);
}

void GraphLaplacian::BuildKDTree() {
  // (re-)build the kd-tree
  kdtree.buildIndex();
}

void GraphLaplacian::FindNeighbors(Eigen::VectorXd const& x, double const r, std::vector<std::pair<std::size_t, double> >& neighbors) const {
  // make sure the state size matches
  assert(x.size()==cloud.StateDim());

  // need to use the squared radius
  const double r2 = r*r;

  // unsorted radius set
  nanoflann::RadiusResultSet<double, std::size_t> resultSet(r2, neighbors); // use squared radius because kd tree metric is the squared euclidean distance

  // find the nearest neighbors---neighbors in a specified radius
  kdtree.findNeighbors(resultSet, x.data(), nanoflann::SearchParams());
}

double GraphLaplacian::FindNeighbors(Eigen::VectorXd const& x, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors) const {
  // make sure the state size matches
  assert(x.size()==cloud.StateDim());

  // make sure we have enough points
  assert(k>0);
  assert(k<=NumSamples());

  // find the nearest neighbors
  std::vector<std::size_t> indices(k);
  std::vector<double> squaredDists(k);
  nanoflann::KNNResultSet<double, std::size_t> resultSet(k);
  resultSet.init(&indices[0], &squaredDists[0]);
  kdtree.findNeighbors(resultSet, x.data(), nanoflann::SearchParams());

  assert(indices.size()==k);
  assert(squaredDists.size()==k);

  // store the indices/distances in the output vector
  neighbors.reserve(indices.size());
  std::transform(indices.begin(), indices.end(), squaredDists.begin(), std::back_inserter(neighbors),
               [](std::size_t a, double b) { return std::make_pair(a, b); });

  return resultSet.worstDist();
}

double GraphLaplacian::EvaluateKernel(Eigen::VectorXd const& x, double const h2, std::vector<std::pair<std::size_t, double> >& neighbors) const {
  assert(kernel);

  // loop through the nearest neighbors
  double sum = 0.0;
  for( auto& neigh : neighbors ) {
    // compute the kernel between the given point and its neighbor
    neigh.second = kernel->EvaluateCompactKernel(neigh.second/h2);
    sum += neigh.second;
  }

  return sum;
}

void GraphLaplacian::ConstructHeatMatrix() {
  // build the kd-tree based on the samples
  BuildKDTree();

  // loop through the samples
  for( std::size_t i=0; i<NumSamples(); ++i ) {
    // get the state for this sample
    const Eigen::Ref<const Eigen::VectorXd> x = Point(i);

    std::cout << x.transpose() << std::endl;
  }

  std::cout << "Construct heat matrix" << std::endl;
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

Eigen::Ref<const Eigen::VectorXd> GraphLaplacian::PointCloud::Point(std::size_t const i) const {
  assert(samples);
  assert(i<samples->size());
  return samples->at(i)->state[0];
}

std::size_t GraphLaplacian::PointCloud::StateDim() const {
  assert(samples);
  assert(samples->size()>0);
  return samples->at(0)->state[0].size();
}
