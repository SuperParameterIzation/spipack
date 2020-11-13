#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

#include <MUQ/Utilities/HDF5/HDF5File.h>

#include "spipack/Tools/Kernels/CompactKernel.hpp"

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

SampleRepresentation::SampleRepresentation(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
samples(std::make_shared<NearestNeighbors>(rv, options["NearestNeighbors"])),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
truncateKernelMatrix(options["TruncateKernelMatrix"].as<bool>(defaults.truncateKernelMatrix)),
numThreads(options["NumThreads"].as<std::size_t>(samples->NumThreads()))
{
  Initialize(options);
}

SampleRepresentation::SampleRepresentation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
samples(std::make_shared<NearestNeighbors>(samples, options["NearestNeighbors"])),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
truncateKernelMatrix(options["TruncateKernelMatrix"].as<bool>(defaults.truncateKernelMatrix)),
numThreads(options["NumThreads"].as<std::size_t>(this->samples->NumThreads()))
{
  Initialize(options);
}

SampleRepresentation::SampleRepresentation(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
samples(samples),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
truncateKernelMatrix(options["TruncateKernelMatrix"].as<bool>(defaults.truncateKernelMatrix)),
numThreads(options["NumThreads"].as<std::size_t>(this->samples->NumThreads()))
{
  Initialize(options);
}

void SampleRepresentation::Initialize(YAML::Node const& options) {
  // create the kernel
  kernel = IsotropicKernel::Construct(options["KernelOptions"]);
  assert(kernel);
  auto compactKernelPtr = std::dynamic_pointer_cast<CompactKernel>(kernel);
  compactKernel = (compactKernelPtr!=nullptr);
  truncationTol = (compactKernel? 1.0 : options["TruncationTolerance"].as<double>(defaults.TruncationTolerance(compactKernel)));
  assert(truncationTol>0.0);
}

std::size_t SampleRepresentation::NumSamples() const { return samples->NumSamples(); }

std::size_t SampleRepresentation::NumNearestNeighbors() const { return numNearestNeighbors; }

Eigen::Ref<Eigen::VectorXd const> SampleRepresentation::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return samples->Point(i);
}

std::shared_ptr<const IsotropicKernel> SampleRepresentation::Kernel() const { return kernel; }

void SampleRepresentation::BuildKDTrees() const { samples->BuildKDTrees(); }

Eigen::VectorXd SampleRepresentation::SquaredBandwidth() const { return samples->SquaredBandwidth(numNearestNeighbors); }

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat) const {
  assert(eps>0.0);

  // compute the squared bandwidth
  const Eigen::VectorXd squaredBandwidth = samples->SquaredBandwidth(numNearestNeighbors);

  // if we are truncating but for some reason want a dense matrix
  if( truncateKernelMatrix ) {
    Eigen::SparseMatrix<double> kmatSparse;
    const Eigen::VectorXd rowsum = KernelMatrix(eps, squaredBandwidth.array().sqrt(), kmatSparse);
    kmat = Eigen::MatrixXd(kmatSparse);
    return rowsum;
  }

  // compute the kernel matrix using the bandwidth
  return KernelMatrix(eps, squaredBandwidth.array().sqrt(), kmat);
}

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, Eigen::Ref<Eigen::MatrixXd> kmat) const {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(rvec.size()==n);

  // check the matrix size
  assert(kmat.rows()==n); assert(kmat.cols()==n);

  // loop through the rows of the kernel matrix
  #pragma omp parallel for num_threads(numThreads)
  for( std::size_t i=0; i<n; ++i ) {
    // compute the bandwith parameter for each pair of points
    Eigen::VectorXd theta = Eigen::VectorXd::Constant(n-i, eps*rvec(i));
    theta = theta.array()*rvec.tail(n-i).array();

    // loop through the columns of the kernel matrix
    for( std::size_t j=i; j<n; ++j ) {
      // evaluate the kernel and insert into the matrix
      const Eigen::VectorXd diff = Point(i)-Point(j);
      kmat(i, j) = kernel->EvaluateIsotropicKernel(diff.dot(diff)/theta(j-i));
      if( j!=i ) { kmat(j, i) = kmat(i, j); }
    }
  }

  // the sum of each row in the kernel matrix
  return kmat.rowwise().sum();
}

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat) const {
  assert(eps>0.0);

  // compute the squared bandwidth
  const Eigen::VectorXd squaredBandwidth = samples->SquaredBandwidth(numNearestNeighbors);

  // compute the kernel matrix using the bandwidth
  if( truncateKernelMatrix ) { return KernelMatrix(eps, squaredBandwidth.array().sqrt(), kmat); }

  // compute a dense  kernel matrix using the bandwidth
  Eigen::MatrixXd kmatDense(NumSamples(), NumSamples());
  const Eigen::VectorXd rowsum = KernelMatrix(eps, squaredBandwidth.array().sqrt(), kmatDense);

  // convert to sparse form
  kmat = kmatDense.sparseView();
  return rowsum;
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
      const double maxband = truncationTol*theta.maxCoeff();
      samples->FindNeighbors(Point(i), maxband, neighbors, i);

      // loop through the neighbors
      for( const auto& neigh : neighbors ) {
        // we have ignored the first i
        assert(neigh.first>=i);

        // skip if we are outside of the kernel's support
        const double para = neigh.second/theta(neigh.first-i);
        if( para>truncationTol ) { continue; }

        // evaluate the kernel
        const double kern = kernel->EvaluateIsotropicKernel(para);

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

void SampleRepresentation::WriteToFile(std::string const& filename, std::string const& dataset) const {
  // create an hdf5 file
  auto file = std::make_shared<HDF5File>(filename);

  // output the collection to file
  const std::string dataset_ = (dataset.at(0)=='/'? dataset : "/"+dataset);
  samples->Samples()->WriteToFile(filename, dataset_);

  file->Close();
}

double SampleRepresentation::DefaultParameters::TruncationTolerance(bool const compact) { return (compact? 1.0 : -std::log(5.0e-2)); }
