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
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim)),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
numThreads(options["NumThreads"].as<std::size_t>(samples->NumThreads()))
{
  Initialize(options);
}

SampleRepresentation::SampleRepresentation(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
samples(std::make_shared<NearestNeighbors>(samples, options["NearestNeighbors"])),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
truncateKernelMatrix(options["TruncateKernelMatrix"].as<bool>(defaults.truncateKernelMatrix)),
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim)),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
numThreads(options["NumThreads"].as<std::size_t>(this->samples->NumThreads()))
{
  Initialize(options);
}

SampleRepresentation::SampleRepresentation(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
samples(samples),
numNearestNeighbors(options["NumNearestNeighbors"].as<std::size_t>(defaults.numNearestNeighbors)),
truncateKernelMatrix(options["TruncateKernelMatrix"].as<bool>(defaults.truncateKernelMatrix)),
manifoldDim(options["ManifoldDimension"].as<double>(defaults.manifoldDim)),
bandwidthPara(options["BandwidthParameter"].as<double>(defaults.bandwidthPara)),
numThreads(options["NumThreads"].as<std::size_t>(this->samples->NumThreads()))
{
  Initialize(options);
}

double SampleRepresentation::BandwidthParameter() const { return bandwidthPara; }

void SampleRepresentation::Initialize(YAML::Node const& options) {
  // create the kernel
  kernel = IsotropicKernel::Construct(options["KernelOptions"]);
  assert(kernel);
  auto compactKernelPtr = std::dynamic_pointer_cast<CompactKernel>(kernel);
  compactKernel = (compactKernelPtr!=nullptr);
  truncationTol = (compactKernel? 1.0 : options["TruncationTolerance"].as<double>(defaults.TruncationTolerance(compactKernel)));
  assert(truncationTol>0.0);

  const YAML::Node& opt = options["Optimization"].as<YAML::Node>(YAML::Node());

  pt.put("Ftol.AbsoluteTolerance", opt["Ftol.AbsoluteTolerance"].as<double>(1.0e-6));
  pt.put("Ftol.RelativeTolerance", opt["Ftol.RelativeTolerance"].as<double>(1.0e-6));
  pt.put("Xtol.AbsoluteTolerance", opt["Xtol.AbsoluteTolerance"].as<double>(1.0e-6));
  pt.put("Xtol.RelativeTolerance", opt["Xtol.RelativeTolerance"].as<double>(1.0e-6));
  pt.put("MaxEvaluations", opt["MaxEvaluations"].as<std::size_t>(1000));
  pt.put("Algorithm", opt["Algorithm"].as<std::string>("COBYLA"));
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

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* options) const {
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
  assert(rvec.rows()==n);

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

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* options) const {
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

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, std::vector<Eigen::Triplet<double> >& entries) const {
    // the number of samples
  const std::size_t n = NumSamples();
  assert(rvec.size()==n);

  // reserve n*log(n) entries (just a guess)
  entries.clear();
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
  return rowsum;
}

Eigen::VectorXd SampleRepresentation::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, Eigen::SparseMatrix<double>& kmat) const {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(rvec.size()==n);

  // compute the entries of the matrix
  std::vector<Eigen::Triplet<double> > entries;
  const Eigen::VectorXd rowsum = KernelMatrix(eps, rvec, entries);

  // resize the matrix
  kmat.resize(n, n);
  kmat.setFromTriplets(entries.begin(), entries.end());

  return rowsum;
}

double SampleRepresentation::KernelDerivativeAverage(double const eps, Eigen::VectorXd const& rvec) const {
  double sm;

  #pragma omp parallel num_threads(numThreads)
  #pragma omp single nowait
  sm = RecursiveKernelDerivativeAverage(0, NumSamples(), eps, rvec, true);

  return sm;
}

double SampleRepresentation::KernelSecondDerivativeAverage(double const eps, Eigen::VectorXd const& rvec) const {
  double sm;

  #pragma omp parallel num_threads(numThreads)
  #pragma omp single nowait
  sm = RecursiveKernelDerivativeAverage(0, NumSamples(), eps, rvec, false);

  return sm;
}

double SampleRepresentation::RecursiveKernelDerivativeRowAverage(std::size_t const row, std::size_t const coli, std::size_t const colj, double const eps, Eigen::VectorXd const& theta, bool const first) const {
  const std::size_t cols = colj-coli;
  const bool includesRow = (coli<=row)&&(row<=colj);

  if( cols<25 ) {
    double sm = 0.0;
    for( std::size_t i=coli; i<colj; ++i ) {
      const Eigen::VectorXd diff = Point(row)-Point(i);
      const double para = diff.dot(diff)/theta(i-row);
      sm -= (i==row? 1.0 : 2.0)*( first? kernel->IsotropicKernelDerivative(para) : -(kernel->IsotropicKernelSecondDerivative(para)*para + 2.0*kernel->IsotropicKernelDerivative(para))/eps )*para/eps;
    }

    return sm/(2*cols-includesRow);
  }

  const std::size_t n = coli+colj;
  const std::size_t half = n/2;
  const bool rowInLeft = includesRow&&(row<half);

  const std::size_t w1 = 2*(half-coli)-rowInLeft;
  const std::size_t w2 = 2*(colj-half)-(!rowInLeft&&includesRow);

  const double x = RecursiveKernelDerivativeRowAverage(row, coli, half, eps, theta, first);
  const double y = RecursiveKernelDerivativeRowAverage(row, half, colj, eps, theta, first);
  return (w1*x + w2*y)/(w1+w2);
}

double SampleRepresentation::KernelDerivativeAverage(std::size_t const rowi, std::size_t const rowj, double const eps, Eigen::VectorXd const& rvec, bool const first) const {
  const std::size_t n = NumSamples();

  // loop through each row
  double sm = 0.0;
  const std::size_t totWeight = (rowj-rowi)*(2*n-rowj-rowi);
  for( std::size_t i=rowi; i<rowj; ++i ) {
    // compute the bandwith parameter for each pair of points
    Eigen::VectorXd theta = Eigen::VectorXd::Constant(n-i, eps*rvec(i));
    theta = theta.array()*rvec.tail(n-i).array();

    if( truncateKernelMatrix ) {
      // find the neighbors within the max bandwidth (ignoring the first i samples)
      std::vector<std::pair<std::size_t, double> > neighbors;
      const double maxband = truncationTol*theta.maxCoeff();
      samples->FindNeighbors(Point(i), maxband, neighbors, i);

      // loop through the neighbors
      double sminner = 0.0;
      for( const auto& neigh : neighbors ) {
        // we have ignored the first i
        assert(neigh.first>=i);

        // skip if we are outside of the kernel's support
        const double para = neigh.second/theta(neigh.first-i);
        if( para>truncationTol ) { continue; }

        sminner -= (i==neigh.first? 1.0 : 2.0)*( first? kernel->IsotropicKernelDerivative(para) : -(kernel->IsotropicKernelSecondDerivative(para)*para + 2.0*kernel->IsotropicKernelDerivative(para))/eps )*para/eps;
      }

      sm += sminner/totWeight;
    } else {
      const double weight = (2*(n-i)-1)/(double)totWeight;
      sm += weight*RecursiveKernelDerivativeRowAverage(i, i, NumSamples(), eps, theta, first);
    }
  }

  return sm;
}

double SampleRepresentation::RecursiveKernelDerivativeAverage(std::size_t const rowi, std::size_t const rowj, double const eps, Eigen::VectorXd const& rvec, bool const first) const {
  assert(rowj>rowi);
  const std::size_t diff = rowj-rowi;
  // if we have to sum over at most 5 rows, use a serial outer loop
  if( diff<5 ) { return KernelDerivativeAverage(rowi, rowj, eps, rvec, first); }

  double x, y;
  const std::size_t tot = rowi+rowj;
  const std::size_t half = tot/2;
  const std::size_t n = NumSamples();

  const std::size_t w1 = (half-rowi)*(2*n-half-rowi);
  const std::size_t w2 = (rowj-half)*(2*n-rowj-half);

  #pragma omp task shared(x)
  x = RecursiveKernelDerivativeAverage(rowi, half, eps, rvec, first);
  #pragma omp task shared(y)
  y = RecursiveKernelDerivativeAverage(half, rowj, eps, rvec, first);
  #pragma omp taskwait
  x = (w1*x+w2*y)/(w1+w2);

  return x;
}

void SampleRepresentation::WriteToFile(std::string const& filename, std::string const& dataset) const {
  // create an hdf5 file
  auto file = std::make_shared<HDF5File>(filename);

  // output the collection to file
  const std::string dataset_ = (dataset.at(0)=='/'? dataset : "/"+dataset);
  samples->Samples()->WriteToFile(filename, dataset_);

  file->Close();
}

double SampleRepresentation::ManifoldDimension() const { return manifoldDim; }

double SampleRepresentation::DefaultParameters::TruncationTolerance(bool const compact) { return (compact? 1.0 : -std::log(5.0e-2)); }
