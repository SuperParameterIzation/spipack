#include "spipack/NumericalSolvers/GraphLaplacian/GraphLaplacian.hpp"

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <MUQ/Utilities/HDF5/HDF5File.h>

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

double GraphLaplacian::DefaultParameters::SquaredBandwidth(YAML::Node const& options) {
  const double bandwidth = options["Bandwidth"].as<double>(defaults.bandwidth);
  return bandwidth*bandwidth;
}

GraphLaplacian::GraphLaplacian(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
  cloud(SampleRandomVariable(rv, options["NumSamples"].as<std::size_t>())),
  kdtree(cloud.StateDim(), cloud, nanoflann::KDTreeSingleIndexAdaptorParams(options["MaxLeaf"].as<std::size_t>(defaults.maxLeaf))),
  bandwidth2(DefaultParameters::SquaredBandwidth(options)),
  eigensolverTol(options["EigensolverTol"].as<double>(defaults.eigensolverTol)),
  eigensolverMaxIt(options["EigensolverMaxIt"].as<std::size_t>(defaults.eigensolverMaxIt))
{
  Initialize(options);
}

GraphLaplacian::GraphLaplacian(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options) :
  cloud(samples),
  kdtree(cloud.StateDim(), cloud, nanoflann::KDTreeSingleIndexAdaptorParams(options["MaxLeaf"].as<std::size_t>(defaults.maxLeaf))),
  bandwidth2(DefaultParameters::SquaredBandwidth(options)),
  eigensolverTol(options["EigensolverTol"].as<double>(defaults.eigensolverTol)),
  eigensolverMaxIt(options["EigensolverMaxIt"].as<std::size_t>(defaults.eigensolverMaxIt))
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

double GraphLaplacian::SquaredBandwidth() const { return bandwidth2; }

std::size_t GraphLaplacian::NumSamples() const { return cloud.kdtree_get_point_count(); }

Eigen::Ref<const Eigen::VectorXd> GraphLaplacian::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return cloud.Point(i);
}

void GraphLaplacian::BuildKDTree() {
  // (re-)build the kd-tree
  kdtree.buildIndex();
}

void GraphLaplacian::FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& x, double const r2, std::vector<std::pair<std::size_t, double> >& neighbors) const {
  // make sure the state size matches
  assert(x.size()==cloud.StateDim());

  // unsorted radius set
  nanoflann::RadiusResultSet<double, std::size_t> resultSet(r2, neighbors); // use squared radius because kd tree metric is the squared euclidean distance

  // find the nearest neighbors---neighbors in a specified radius
  kdtree.findNeighbors(resultSet, x.data(), nanoflann::SearchParams());
}

double GraphLaplacian::FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& x, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors) const {
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

double GraphLaplacian::EvaluateKernel(Eigen::Ref<const Eigen::VectorXd> const& x, double const h2, std::vector<std::pair<std::size_t, double> >& neighbors) const {
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

  // the number of samples
  const std::size_t n = NumSamples();

  // the sum of the kernel functions between sample i and its nearest neighbors
  Eigen::VectorXd kernelsum(n);

  // loop through the samples
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors(n);
  std::size_t numentries = 0; // used to count the number of nonzero entries
  for( std::size_t i=0; i<n; ++i ) {
    // get the state for this sample
    const Eigen::Ref<const Eigen::VectorXd> x = Point(i);

    // find the nearest neighbors
    FindNeighbors(x, bandwidth2, neighbors[i]);
    double h2 = bandwidth2;
    // make sure we have at least 3 neighbors
    if( neighbors[i].size()<3 ) {
      neighbors[i].clear();
      h2 = FindNeighbors(x, (std::size_t)3, neighbors[i]);
    }

    numentries += neighbors[i].size();

    // evaluate the kernel function at the neighbors
    kernelsum(i) = EvaluateKernel(x, h2, neighbors[i]);
  }

  // construct the nonzero entries of the heat matrix
  ConstructHeatMatrix(kernelsum, neighbors, numentries);
}

void GraphLaplacian::ConstructHeatMatrix(Eigen::Ref<const Eigen::VectorXd> const& kernelsum, std::vector<std::vector<std::pair<std::size_t, double> > >& neighbors, std::size_t const numentries) {
  // create the nonzero entries
  std::vector<Eigen::Triplet<double> > entries;
  entries.reserve(numentries);

  // the number of samples
  const std::size_t n = NumSamples();

  // loop through each element
  for( std::size_t i=0; i<n; ++i ) {
    // loop through the neighbors and comptue the corrected kernel
    double sum = 0.0;
    for( auto& it : neighbors[i] ) {
      it.second /= std::sqrt(kernelsum(i)*kernelsum(it.first));
      sum += it.second;
    }

    // loop through the neighbors and add the normalized kernel to the matrix
    for( const auto& it : neighbors[i] ) {
      entries.push_back(Eigen::Triplet<double>(i, it.first, it.second/sum));
    }
  }

  // create the heat matrix
  heatMatrix.setFromTriplets(entries.begin(), entries.end());
}

Eigen::VectorXd GraphLaplacian::HeatMatrixEigenvalues(const size_t neig) const {
  return ComputeLargestSparseEigenvalues(neig, heatMatrix);
}

Eigen::VectorXd GraphLaplacian::ComputeSparseEigenvalues(std::size_t const neig, Eigen::SparseMatrix<double> const& mat, bool const computeLargest) const {
  if( computeLargest ) { return ComputeLargestSparseEigenvalues(neig, mat); }
  return ComputeSmallestSparseEigenvalues(neig, mat);
}

Eigen::VectorXd GraphLaplacian::ComputeLargestSparseEigenvalues(std::size_t const neig, Eigen::SparseMatrix<double> const& mat) const {
  // wrapper for space mat-vecs
  Spectra::SparseGenMatProd<double> matvec(heatMatrix);

  // construct eigen solver object, requesting the largest neig eigenvalues
  Spectra::GenEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolver(&matvec, neig, 10*neig);

  // initialize and compute
  eigsolver.init();
  const int ncomputed = eigsolver.compute(eigensolverMaxIt, eigensolverTol);
  assert(ncomputed==neig);

  std::cout << "largest: " << eigsolver.eigenvalues() << std::endl;

  // get the results
  return eigsolver.eigenvalues().real();
}

Eigen::VectorXd GraphLaplacian::ComputeSmallestSparseEigenvalues(std::size_t const neig, Eigen::SparseMatrix<double> const& mat) const {
  // wrapper for space mat-vecs
  Spectra::SparseGenMatProd<double> matvec(heatMatrix);

  // construct eigen solver object, requesting the largest neig eigenvalues
  Spectra::GenEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolver(&matvec, neig, 10*neig);

  // initialize and compute
  eigsolver.init();
  const int ncomputed = eigsolver.compute(eigensolverMaxIt, eigensolverTol);
  std::cout << "smallest: " << eigsolver.eigenvalues().transpose() << std::endl;
  assert(ncomputed==neig);

  // get the results
  return eigsolver.eigenvalues().real();
}

const Eigen::Ref<const Eigen::SparseMatrix<double> > GraphLaplacian::HeatMatrix() const { return heatMatrix; }

void GraphLaplacian::SolveWeightedPoisson(Eigen::Ref<Eigen::VectorXd> vec) {
  // construct the heat matrix
  ConstructHeatMatrix();

  std::cout << "num heat matrix entries: " << heatMatrix.nonZeros() << std::endl;

  // the number of samples
  const std::size_t n = NumSamples();
  assert(n==heatMatrix.rows());
  assert(n==heatMatrix.cols());

  // initalize the discrete Laplacian as an identity
  Eigen::SparseMatrix<double> laplace(n, n);
  laplace.setIdentity();

  std::cout << "laplace matrix entries (identity): " << laplace.nonZeros() << std::endl;

  laplace -= heatMatrix;
  laplace /= bandwidth2;

  std::cout << "laplace matrix entries (no bc): " << laplace.nonZeros() << std::endl;

  //std::cout << "smallest eig values (no bc): " << ComputeSmallestSparseEigenvalues(1, laplace).transpose() << std::endl;

  //std::cout << "largest eig values (no bc): " << ComputeLargestSparseEigenvalues(1, laplace).transpose() << std::endl;

  /*for( std::size_t j=0; j<2; ++j ) {
    for( std::size_t i=0; i<n; ++i ) {
      std::cout << laplace.coeff(j, i) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/

  //if( laplace.isCompressed() ) { laplace.uncompress(); }

  //laplace.prune([](int i, int j, float) { return (i!=0 && i!=1); });
  laplace.prune([](int i, int j, float) { return i!=0; });
  for( std::size_t i=0; i<n; ++i ) {
    laplace.coeffRef(0, i) = 1.0;
    //laplace.coeffRef(1, i) = 1.0;
  }
  //laplace.coeffRef(0, 0) = 1.0;
  //laplace.coeffRef(1, 1) = 1.0;
  vec(0) = 0.0;
  //vec(1) = 0.0;

  /*for( std::size_t j=0; j<2; ++j ) {
    for( std::size_t i=0; i<n; ++i ) {
      std::cout << laplace.coeff(j, i) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;*/

  std::cout << "laplace matrix entries (with bc): " << laplace.nonZeros() << std::endl;


  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  // Compute the ordering permutation vector from the structural pattern of A
solver.analyzePattern(laplace);
// Compute the numerical factorization
solver.factorize(laplace);
  assert(solver.info()==Eigen::Success);
  auto vec0 = solver.solve(vec);

  std::cout << "sollved!" << std::endl;

  //std::cout << vec0.transpose() << std::endl;
}

double GraphLaplacian::EigensolverTolerance() const { return eigensolverTol; }

void GraphLaplacian::WriteToFile(std::string const& filename, std::string const& dataset) const {
  // create an hdf5 file
  auto file = std::make_shared<HDF5File>(filename);

  // output the collection to file
  const std::string dataset_ = (dataset.at(0)=='/'? dataset : "/"+dataset);
  cloud.samples->WriteToFile(filename, dataset);

  file->Close();
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
