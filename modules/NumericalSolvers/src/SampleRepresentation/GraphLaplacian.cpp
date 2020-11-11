#include "spipack/NumericalSolvers/SampleRepresentation/GraphLaplacian.hpp"

#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include "spipack/Tools/SortVector.hpp"

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

GraphLaplacian::GraphLaplacian(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
bandwidthRange(std::pair<double, double>(options["BandwidthRange.Min"].as<double>(defaults.bandwidthRange.first),options["BandwidthRange.Max"].as<double>(defaults.bandwidthRange.second))),
numBandwidthSteps(options["NumBandwidthSteps"].as<std::size_t>(defaults.numBandwidthSteps)),
tuneBandwidthParameter(options["TuneBandwidth"].as<bool>(defaults.tuneBandwidthParameter)),
eigensolverTol(options["EigensolverTol"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIt"].as<std::size_t>(defaults.eigensolverMaxIt))
{
  Initialize(options);
}

GraphLaplacian::GraphLaplacian(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
bandwidthRange(std::pair<double, double>(options["BandwidthRange.Min"].as<double>(defaults.bandwidthRange.first),options["BandwidthRange.Max"].as<double>(defaults.bandwidthRange.second))),
numBandwidthSteps(options["NumBandwidthSteps"].as<std::size_t>(defaults.numBandwidthSteps)),
tuneBandwidthParameter(options["TuneBandwidth"].as<bool>(defaults.tuneBandwidthParameter)),
eigensolverTol(options["EigensolverTol"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIt"].as<std::size_t>(defaults.eigensolverMaxIt))
{
  Initialize(options);
}

void GraphLaplacian::Initialize(YAML::Node const& options) {
  // initialize the sparse heat matrix
  heatMatrix.resize(NumSamples(), NumSamples());
}

bool GraphLaplacian::TuneBandwidthParameter() const { return tuneBandwidthParameter; }

void GraphLaplacian::BuildKDTrees() {
  // (re-)build the kd-tree
  samples.BuildKDTrees();
}

Eigen::VectorXd GraphLaplacian::SquaredBandwidth(std::vector<std::vector<std::pair<std::size_t, double> > >& neighbors) const {
  // the number of samples
  const std::size_t n = NumSamples();

  // loop through each sample and compute the squared bandwidth
  Eigen::VectorXd bandwidth(n);
  neighbors.resize(n);
  for( std::size_t i=0; i<n; ++i ) {
    // find the nearest neighbors for each sample
    bandwidth(i) = samples.FindNeighbors(Point(i), numNearestNeighbors, neighbors[i]);
  }

  return bandwidth;
}

Eigen::VectorXd GraphLaplacian::SquaredBandwidth() const {
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  return SquaredBandwidth(neighbors);
}

Eigen::VectorXd GraphLaplacian::KernelMatrix(double const bandwidthPara, Eigen::Ref<Eigen::VectorXd const> const& bandwidth, std::vector<std::vector<std::pair<std::size_t, double> > > const& neighbors, Eigen::SparseMatrix<double>& kernmat) const {
  assert(kernmat.cols()==bandwidth.size());
  assert(kernmat.rows()==bandwidth.size());

  // reserve n*log(n) entries (just a guess)
  std::vector<Eigen::Triplet<double> > entries;
  entries.reserve(std::size_t(NumSamples()*std::log((double)NumSamples())));

  Eigen::VectorXd rowsum = Eigen::VectorXd::Zero(bandwidth.size());

  for( std::size_t i=0; i<bandwidth.size(); ++i ) {
    // compute the bandwith parameter for each pair of points
    Eigen::VectorXd theta = Eigen::VectorXd::Constant(bandwidth.size()-i, bandwidthPara*bandwidth(i));
    theta = theta.array()*bandwidth.tail(bandwidth.size()-i).array();

    for( const auto& neigh : neighbors[i] ) {
      if( neigh.first<i ) { continue; }

      const double para = neigh.second/theta(neigh.first-i);
      if( para>1.0 ) { continue; }

      // evaluate the kernel
      const double kern = kernel->EvaluateCompactKernel(para);

      // insert into the matrix
      entries.push_back(Eigen::Triplet<double>(i, neigh.first, kern));
      rowsum(i) += kern;
      if( neigh.first!=i ) {
        entries.push_back(Eigen::Triplet<double>(neigh.first, i, kern));
        rowsum(i) += kern;
      }
    }
  }

  kernmat.setFromTriplets(entries.begin(), entries.end());
  return rowsum;
}

Eigen::VectorXd GraphLaplacian::KernelMatrix(double const bandwidthPara, Eigen::Ref<Eigen::VectorXd const> const& bandwidth, Eigen::SparseMatrix<double>& kernmat) const {
  assert(kernmat.cols()==bandwidth.size());
  assert(kernmat.rows()==bandwidth.size());

  // reserve n*log(n) entries (just a guess)
  std::vector<Eigen::Triplet<double> > entries;
  entries.reserve(std::size_t(NumSamples()*std::log((double)NumSamples())));

  Eigen::VectorXd rowsum = Eigen::VectorXd::Zero(bandwidth.size());

  for( std::size_t i=0; i<bandwidth.size(); ++i ) {
    // compute the bandwith parameter for each pair of points
    Eigen::VectorXd theta = Eigen::VectorXd::Constant(bandwidth.size()-i, bandwidthPara*bandwidth(i));
    theta = theta.array()*bandwidth.tail(bandwidth.size()-i).array();

    // the max bandwidth
    std::vector<std::pair<std::size_t, double> > neighbors;
    const double maxband = theta.maxCoeff();
    samples.FindNeighbors(Point(i), maxband, neighbors, i);

    for( const auto& neigh : neighbors ) {
      assert(neigh.first>=i);

      const double para = neigh.second/theta(neigh.first-i);
      if( para>1.0 ) { continue; }

      // evaluate the kernel
      const double kern = kernel->EvaluateCompactKernel(para);

      // insert into the matrix
      entries.push_back(Eigen::Triplet<double>(i, neigh.first, kern));
      rowsum(i) += kern;
      if( neigh.first!=i ) {
        entries.push_back(Eigen::Triplet<double>(neigh.first, i, kern));
        rowsum(i) += kern;
      }
    }
  }

  kernmat.setFromTriplets(entries.begin(), entries.end());
  return rowsum;
}

Eigen::VectorXd GraphLaplacian::DensityEstimation() const {
  // compute the bandwidth
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  const Eigen::VectorXd squaredBandwidth = SquaredBandwidth(neighbors);

  if( tuneBandwidthParameter ) {
    std::cout << "tune" << std::endl;
  } else {
    std::cout << "do not tune" << std::endl;
  }

  return Eigen::VectorXd();
}

Eigen::Matrix<double, Eigen::Dynamic, 2> GraphLaplacian::TuneKernelBandwidth(Eigen::Ref<Eigen::VectorXd const> const& bandwidth, std::vector<std::vector<std::pair<std::size_t, double> > > const& neighbors, std::pair<double, Eigen::SparseMatrix<double> >& optimal) const {
  // get the candidate bandwidth parameters
  const Eigen::VectorXd bandwidthPara = BandwidthParameterCandidates();

  // the number of samples
  const std::size_t n = NumSamples();

  // the potential kernel matrices
  auto kernmatBest = std::make_shared<Eigen::SparseMatrix<double> >(n, n);

  // loop through the possible bandwith parameters
  Eigen::Matrix<double, Eigen::Dynamic, 2> sigmaprime(bandwidthPara.size()-1, 2);
  double logsig, logsigprev;
  std::size_t ind;
  double maxsigprime = 0.0;
  for( std::size_t l=0; l<bandwidthPara.size(); ++l ) {
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "eps: " << bandwidthPara(l) << std::endl;
    auto kernmat = std::make_shared<Eigen::SparseMatrix<double> >(n, n);
    const Eigen::VectorXd rowsum = KernelMatrix(bandwidthPara(l), bandwidth, neighbors, *kernmat);
    std::cout << "computed kern mat: " << kernmat->nonZeros() << " of " << n*n << std::endl;

    logsigprev = logsig;
    logsig = std::log(rowsum.sum()/((double)n*n));

    std::cout << "comptued sum" << std::endl;
    if( l>0 ) {
      sigmaprime(l-1, 0) = bandwidthPara(l-1);
      sigmaprime(l-1, 1) = (logsig-logsigprev)/(std::log(bandwidthPara(l))-std::log(bandwidthPara(l-1)));
      if( sigmaprime(l-1, 1)>maxsigprime ) {
        ind = l-1;
        maxsigprime = sigmaprime(l-1, 1);
        kernmatBest = kernmat;
      }
    }
  }

  // set the output
  optimal.first = bandwidthPara(ind);
  optimal.second = *kernmatBest;

  return sigmaprime;
}

Eigen::Matrix<double, Eigen::Dynamic, 2> GraphLaplacian::TuneKernelBandwidth(Eigen::Ref<Eigen::VectorXd const> const& bandwidth, std::pair<double, Eigen::SparseMatrix<double> >& optimal) const {
  // get the candidate bandwidth parameters
  const Eigen::VectorXd bandwidthPara = BandwidthParameterCandidates();

  // the number of samples
  const std::size_t n = NumSamples();

  // the potential kernel matrices
  auto kernmatBest = std::make_shared<Eigen::SparseMatrix<double> >(n, n);

  // loop through the possible bandwith parameters
  Eigen::Matrix<double, Eigen::Dynamic, 2> sigmaprime(bandwidthPara.size()-1, 2);
  double logsig, logsigprev;
  std::size_t ind;
  double maxsigprime = 0.0;
  for( std::size_t l=0; l<bandwidthPara.size(); ++l ) {
    std::cout << std::endl << std::endl << std::endl;
    std::cout << "eps: " << bandwidthPara(l) << std::endl;
    auto kernmat = std::make_shared<Eigen::SparseMatrix<double> >(n, n);
    const Eigen::VectorXd rowsum = KernelMatrix(bandwidthPara(l), bandwidth, *kernmat);
    std::cout << "computed kern mat: " << kernmat->nonZeros() << " of " << n*n << std::endl;

    logsigprev = logsig;
    logsig = std::log(rowsum.sum()/((double)n*n));

    std::cout << "comptued sum" << std::endl;
    if( l>0 ) {
      sigmaprime(l-1, 0) = bandwidthPara(l-1);
      sigmaprime(l-1, 1) = (logsig-logsigprev)/(std::log(bandwidthPara(l))-std::log(bandwidthPara(l-1)));
      if( sigmaprime(l-1, 1)>maxsigprime ) {
        ind = l-1;
        maxsigprime = sigmaprime(l-1, 1);
        kernmatBest = kernmat;
      }
    }
  }

  // set the output
  optimal.first = bandwidthPara(ind);
  optimal.second = *kernmatBest;

  return sigmaprime;
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

std::size_t GraphLaplacian::NumBandwidthSteps() const { return numBandwidthSteps; }

std::pair<double, double> GraphLaplacian::BandwidthRange() const { return bandwidthRange; }

Eigen::VectorXd GraphLaplacian::BandwidthParameterCandidates() const {
  // the candidate exponents
  Eigen::VectorXd candidates = Eigen::VectorXd::LinSpaced(numBandwidthSteps+1, bandwidthRange.first, bandwidthRange.second);

  // raise 2 to the exponent
  for( std::size_t i=0; i<candidates.size(); ++i ) { candidates(i) = std::pow(2.0, candidates(i)); }

  return candidates;
}

void GraphLaplacian::ConstructHeatMatrix() {
  /*// build the kd-tree based on the samples
  BuildKDTrees();

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
    samples.FindNeighbors(x, bandwidth2, neighbors[i]);
    double h2 = bandwidth2;
    // make sure we have at least 3 neighbors
    if( neighbors[i].size()<3 ) {
      neighbors[i].clear();
      h2 = samples.FindNeighbors(x, (std::size_t)3, neighbors[i]);
    }

    numentries += neighbors[i].size();

    // evaluate the kernel function at the neighbors
    kernelsum(i) = EvaluateKernel(x, h2, neighbors[i]);
  }

  // construct the nonzero entries of the heat matrix
  ConstructHeatMatrix(kernelsum, neighbors, numentries);*/
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
  //laplace /= bandwidth2;

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
