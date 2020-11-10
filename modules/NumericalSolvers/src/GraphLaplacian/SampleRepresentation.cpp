#include "spipack/NumericalSolvers/GraphLaplacian/SampleRepresentation.hpp"

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

SampleRepresentation::SampleRepresentation(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
samples(rv, options["NearestNeighbors"]),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
numThreads(options["NumThreads"].as<std::size_t>(samples.NumThreads()))
{
  // create the kernel
  kernel = CompactKernel::Construct(options["KernelOptions"]);
  assert(kernel);
}

SampleRepresentation::SampleRepresentation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
samples(samples, options["NearestNeighbors"]),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
numThreads(options["NumThreads"].as<std::size_t>(this->samples.NumThreads()))
{
  // create the kernel
  kernel = CompactKernel::Construct(options["KernelOptions"]);
  assert(kernel);
}

std::size_t SampleRepresentation::NumSamples() const { return samples.NumSamples(); }

std::size_t SampleRepresentation::NumNearestNeighbors() const { return numNearestNeighbors; }

Eigen::Ref<Eigen::VectorXd const> SampleRepresentation::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return samples.Point(i);
}

std::shared_ptr<const CompactKernel> SampleRepresentation::Kernel() const { return kernel; }

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat) const {
  assert(eps>0.0);

  // construct the kd-trees
  samples.BuildKDTrees();

  // compute the squared bandwidth
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  const Eigen::VectorXd squaredBandwidth = samples.SquaredBandwidth(numNearestNeighbors, neighbors);

  // compute the kernel matrix using the bandwidth
  return KernelMatrix(eps, squaredBandwidth.array().sqrt(), kmat);
}

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, Eigen::SparseMatrix<double>& kmat) const {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(rvec.size()==n);

  // resize the matrix
  kmat.resize(n, n);

  // reserve n*log(n) entries (just a guess)
  std::vector<Eigen::Triplet<double> > entries;
  entries.reserve(std::max((std::size_t)1, (std::size_t)(n*std::log((double)n))));

  #pragma omp parallel num_threads(numThreads)
  {
    std::vector<Eigen::Triplet<double> > entries_private;
    entries_private.reserve(std::max((std::size_t)1, (std::size_t)(n/numThreads*std::log((double)n))));

    // loop through the rows of the kernel matrix
    #pragma omp for nowait
    for( std::size_t i=0; i<n; ++i ) {
      // compute the bandwith parameter for each pair of points
      Eigen::VectorXd theta = Eigen::VectorXd::Constant(n-i, eps*rvec(i));
      theta = theta.array()*rvec.tail(n-i).array();

      // find the neighbors within the max bandwidth (ignoring the first i samples)
      std::vector<std::pair<std::size_t, double> > neighbors;
      const double maxband = theta.maxCoeff();
      samples.FindNeighbors(Point(i), maxband, neighbors, i);

      // loop through the neighbors
      for( const auto& neigh : neighbors ) {
        // we have ignored the first i
        assert(neigh.first>=i);

        // skip if we are outside of the kernel's support
        const double para = neigh.second/theta(neigh.first-i);
        if( para>1.0 ) { continue; }

        // evaluate the kernel
        const double kern = kernel->EvaluateCompactKernel(para);

        // insert into the matrix and increment the rowsum
        entries_private.push_back(Eigen::Triplet<double>(i, neigh.first, kern));
        if( neigh.first!=i ) {
          entries_private.push_back(Eigen::Triplet<double>(neigh.first, i, kern));
        }
      }
    }

    #pragma omp critical
    entries.insert(entries.end(), std::make_move_iterator(entries_private.begin()), std::make_move_iterator(entries_private.end()));
  } // end parallel

  // the sum of each row in the kernel matrix
  Eigen::VectorXd rowsum = Eigen::VectorXd::Zero(n);
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }

  kmat.setFromTriplets(entries.begin(), entries.end());
  return rowsum;
}
