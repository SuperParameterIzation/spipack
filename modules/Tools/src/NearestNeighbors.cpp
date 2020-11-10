#include "spipack/Tools/NearestNeighbors.hpp"

#include <omp.h>

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;

NearestNeighbors::NearestNeighbors(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
samples(SampleRandomVariable(rv, options["NumSamples"].as<std::size_t>())),
numThreads(options["NumThreads"].as<std::size_t>(defaults.numThreads))
{
  assert(samples);
  Initialize(options);
}

NearestNeighbors::NearestNeighbors(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
samples(samples),
numThreads(options["NumThreads"].as<std::size_t>(defaults.numThreads))
{
  assert(samples);
  Initialize(options);
}

void NearestNeighbors::Initialize(YAML::Node const& options) {
  // get the stride
  const std::size_t stride = options["Stride"].as<double>(NumSamples()<5? NumSamples() : (size_t)(NumSamples()/5));

  // reserve memory
  clouds.reserve(NumSamples()/stride+1);
  kdtrees.reserve(clouds.size());

  // create a bunch of point clouds
  for( std::size_t lag=0; lag<NumSamples(); lag+=stride ) {
    clouds.emplace_back(samples, lag);
    kdtrees.push_back(
      std::make_shared<NanoflannKDTree>(
        StateDim(),
        *(clouds.end()-1), nanoflann::KDTreeSingleIndexAdaptorParams(options["MaxLeaf"].as<std::size_t>(defaults.maxLeaf))
      ));
  }
}

std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> NearestNeighbors::SampleRandomVariable(std::shared_ptr<RandomVariable> const& rv, std::size_t const n) {
  // add random samples into a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  return samples;
}

Eigen::Ref<Eigen::VectorXd const> NearestNeighbors::Point(std::size_t const i) const {
  assert(samples);
  assert(i<samples->size());
  return samples->at(i)->state[0];
}

std::size_t NearestNeighbors::NumSamples() const {
  assert(samples);
  return samples->size();
}

std::size_t NearestNeighbors::NumThreads() const { return numThreads; }

std::size_t NearestNeighbors::StateDim() const {
  assert(samples);
  assert(samples->size()>0);
  return samples->at(0)->state[0].size();
}

std::shared_ptr<const SampleCollection> NearestNeighbors::Samples() const {
  assert(samples);
  return samples;
}

void NearestNeighbors::BuildKDTrees() const {
  #pragma omp parallel for num_threads(numThreads)
  for( const auto& kdtree : kdtrees ) {
    assert(kdtree);
    // (re-)build the kd-tree
    kdtree->buildIndex();
  }
}

void NearestNeighbors::FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, double const radius2, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag) const {
  // make sure the state size matches
  assert(point.size()==StateDim());

  // unsorted radius set
  nanoflann::RadiusResultSet<double, std::size_t> resultSet(radius2, neighbors); // use squared radius because kd tree metric is the squared euclidean distance

  // choose the kdtree that ignores the first lag samples
  std::size_t ind = clouds.size()-1;
  while( clouds[ind].lag>lag ) { --ind; }

  // find the nearest neighbors---neighbors in a specified radius
  kdtrees[ind]->findNeighbors(resultSet, point.data(), nanoflann::SearchParams());

  // add the lag back to the global indcies
  for( auto& neigh : neighbors ) {
    neigh.first += clouds[ind].lag;
  }

  // remove the ones that should have been ignored
  auto it = std::remove_if(neighbors.begin(), neighbors.end(), [lag](std::pair<std::size_t, double> const& neigh) { return neigh.first<lag; } );
  neighbors.erase(it, neighbors.end());
}

double NearestNeighbors::FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag) const {
  // make sure the state size matches
  assert(point.size()==StateDim());
  assert(k>0);

  // choose the kdtree that ignores the first lag samples
  std::size_t ind = clouds.size()-1;
  while( clouds[ind].lag>lag ) { --ind; }

  // find the nearest neighbors---neighbors in a specified radius
  std::vector<std::size_t> neighInd(k);
  std::vector<double> neighDist(k);
  const std::size_t nfound = kdtrees[ind]->knnSearch(point.data(), k, neighInd.data(), neighDist.data());

  // add the lag back to the global indcies
  neighbors.reserve(k);
  for( std::size_t i=0; i<nfound; ++i ) {
    neighbors.push_back(std::pair<std::size_t, double>(neighInd[i]+clouds[ind].lag, neighDist[i]));
  }

  // remove the ones that should have been ignored
  auto it = std::remove_if(neighbors.begin(), neighbors.end(), [lag](std::pair<std::size_t, double> const& neigh) { return neigh.first<lag; } );
  neighbors.erase(it, neighbors.end());

  double avg = 0.0;
  int sub = 0 ;
  for( const auto& neigh : neighbors ) {
    if( neigh.second<1.0e-12 ) { ++sub; continue; } // don't include itself
    avg += neigh.second;
  }
  return (sub==neighbors.size()? 0.0 : avg/(neighbors.size()-sub));
}

Eigen::VectorXd NearestNeighbors::SquaredBandwidth(std::size_t const k, std::vector<std::vector<std::pair<std::size_t, double> > >& neighbors) const {
  // the number of samples
  const std::size_t n = NumSamples();

  // need to find k+1 neighbors because one will be itself
  const std::size_t kp1 = k+1;

  // loop through each sample and compute the squared bandwidth
  Eigen::VectorXd bandwidth(n);
  neighbors.resize(n);
  #pragma omp parallel for num_threads(numThreads)
  for( std::size_t i=0; i<n; ++i ) {
    // find the nearest neighbors for each sample
    bandwidth(i) = FindNeighbors(Point(i), kp1, neighbors[i]);
  }

  return bandwidth;
}

NearestNeighbors::Cloud::Cloud(std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> const& samples, std::size_t const lag) :
samples(samples),
lag(lag)
{
  assert(samples->size()>lag);
}

std::size_t NearestNeighbors::Cloud::kdtree_get_point_count() const {
  assert(samples);
  return samples->size()-lag;
}

double  NearestNeighbors::Cloud::kdtree_get_pt(std::size_t const p, std::size_t const i) const {
  assert(samples);
  return samples->at(lag+p)->state[0] [i];
}
