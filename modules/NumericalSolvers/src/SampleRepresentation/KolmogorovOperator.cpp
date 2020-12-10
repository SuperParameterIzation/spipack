#include "spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <MUQ/Optimization/NLoptOptimizer.h>

#include "spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<RandomVariable> const& rv, YAML::Node const& options) :
SampleRepresentation(rv, options),
density(std::make_shared<DensityEstimation>(this->samples, DensityOptions(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
exponentPara(options["BandwidthExponent"].as<double>(defaults.exponentPara)),
neig(options["NumEigenvalues"].as<std::size_t>(defaults.neig)),
eigensolverTol(options["EigensolverTolerance"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIterations"].as<std::size_t>(defaults.eigensolverMaxIt))
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<SampleCollection> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(std::make_shared<DensityEstimation>(this->samples, DensityOptions(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
exponentPara(options["BandwidthExponent"].as<double>(defaults.exponentPara)),
neig(options["NumEigenvalues"].as<std::size_t>(defaults.neig)),
eigensolverTol(options["EigensolverTolerance"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIterations"].as<std::size_t>(defaults.eigensolverMaxIt))
{}

KolmogorovOperator::KolmogorovOperator(std::shared_ptr<const NearestNeighbors> const& samples, YAML::Node const& options) :
SampleRepresentation(samples, options),
density(std::make_shared<DensityEstimation>(this->samples, DensityOptions(options))), // default to using the same parameters as the Kolmogorov operator
operatorConstant(options["OperatorParameter"].as<double>(defaults.operatorConstant)),
exponentPara(options["BandwidthExponent"].as<double>(defaults.exponentPara)),
neig(options["NumEigenvalues"].as<std::size_t>(defaults.neig)),
eigensolverTol(options["EigensolverTolerance"].as<double>(defaults.eigensolverTol)),
eigensolverMaxIt(options["EigensolverMaxIterations"].as<std::size_t>(defaults.eigensolverMaxIt))
{}

YAML::Node KolmogorovOperator::DensityOptions(YAML::Node const& options) {
  if( YAML::Node parameter = options["DensityOptions"] ) { return parameter; }
  return options;
}

Eigen::VectorXd KolmogorovOperator::EstimateDensity(bool const tune) const {
  assert(density);

  // estimate the density at each sample
  const Eigen::VectorXd dens = density->Estimate(tune);

  // the tuning (may have) changed the manifold dimension estimation
  manifoldDim = density->ManifoldDimension();

  return dens;
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* tune) {
  assert(kmat.rows()==NumSamples());
  assert(kmat.cols()==NumSamples());

  // estimate the density at each sample
  const Eigen::VectorXd dens = EstimateDensity(*(const bool*)tune);

  // if we are truncating but for some reason want a dense matrix
  if( truncateKernelMatrix ) {
    Eigen::SparseMatrix<double> kmatSparse;
    const Eigen::MatrixXd rowsum = KernelMatrix(eps, dens, kmatSparse);
    kmat = Eigen::MatrixXd(kmatSparse);
    return rowsum;
  }

  return KernelMatrix(eps, dens, kmat);
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::Ref<Eigen::MatrixXd> kmat) {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(dens.rows()==n);

  // compute (and store) the diagonal matrix P
  P = dens.array().pow(exponentPara);

  // compute the unnormalized kernel matrix K
  Eigen::MatrixXd rowsum = SampleRepresentation::KernelMatrix(4.0*eps, P, kmat);

  const double para = VariableBandwidthExponent();

  // compute the normalization for the new kernel matrix
  rowsum = rowsum.array()/P.array().pow(manifoldDim*para);

  // renormalize the kernel matrix
  #pragma omp parallel for num_threads(numThreads)
  for( std::size_t i=0; i<n; ++i ) {
    for( std::size_t j=0; j<n; ++j ) { kmat(i, j) /= rowsum(i)*rowsum(j); }
  }

  return kmat.rowwise().sum();
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* tune) {
  // estimate the density at each sample
  const Eigen::VectorXd dens = EstimateDensity(*(const bool*)tune);

  if( truncateKernelMatrix ) { return KernelMatrix(eps, dens, kmat); }

  // compute a dense kernel matrix using the bandwidth
  Eigen::MatrixXd kmatDense(NumSamples(), NumSamples());
  const Eigen::MatrixXd rowsum = KernelMatrix(eps, dens, kmatDense);

  // convert to sparse form
  kmat = kmatDense.sparseView();
  return rowsum;
}

Eigen::VectorXd KolmogorovOperator::KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::SparseMatrix<double>& kmat) {
  // the number of samples
  const std::size_t n = NumSamples();
  assert(dens.size()==n);

  const double para = VariableBandwidthExponent();

  // compute (and store) the diagonal matrix P
  P = dens.array().pow(exponentPara);

  // compute the unnormalized kernel matrix K
  std::vector<Eigen::Triplet<double> > entries;
  Eigen::MatrixXd rowsum = SampleRepresentation::KernelMatrix(4.0*eps, P, entries);

  // compute the normalization for the new kernel matrix
  rowsum = rowsum.array()/P.array().pow(manifoldDim*para);

  // renormalize the kernel matrix
  #pragma omp parallel for num_threads(numThreads)
  for( auto& entry : entries ) {
    entry = Eigen::Triplet(entry.row(), entry.col(), entry.value()/(rowsum(entry.row())*rowsum(entry.col())));
  }

  // resize the matrix
  kmat.resize(n, n);
  kmat.setFromTriplets(entries.begin(), entries.end());

  // recompute the rowsum
  rowsum = Eigen::VectorXd::Zero(n);
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }
  return rowsum;
}

void KolmogorovOperator::TuneBandwidthParameter(bool const tuneDens) {
  const Eigen::VectorXd rho = EstimateDensity(tuneDens).array().pow(exponentPara);
  manifoldDim = Density()->ManifoldDimension();

  // the cost function and optimizater
  auto cost = std::make_shared<BandwidthCost>(rho, Density());
  auto opt = std::make_shared<NLoptOptimizer>(cost, pt);

  // the initial condition for the optimization is the current parameter value
  std::vector<Eigen::VectorXd> inputs(1);
  assert(bandwidthPara>0.0);
  inputs[0] = Eigen::VectorXd::Constant(1, std::log2(4.0*bandwidthPara));

  // solve the optimization and update the parameters
  std::pair<Eigen::VectorXd, double> soln = opt->Solve(inputs);
  bandwidthPara = std::pow(2, soln.first(0))/4.0;
}

double KolmogorovOperator::VariableBandwidthExponent() const { return 1.0+0.5*manifoldDim*exponentPara+exponentPara-0.5*operatorConstant; }

double KolmogorovOperator::ExponentParameter() const { return exponentPara; }

std::shared_ptr<DensityEstimation> KolmogorovOperator::Density() const { return density; }

void KolmogorovOperator::ComputeEigendecomposition(Eigen::Ref<Eigen::VectorXd> S, Eigen::Ref<Eigen::VectorXd> Sinv, Eigen::Ref<Eigen::VectorXd> eigenvalues, Eigen::Ref<Eigen::MatrixXd> eigenvectors) {
  // the number of samples
  const std::size_t n = NumSamples();

  assert(eigenvalues.size()==neig);
  assert(eigenvectors.rows()==n);
  assert(eigenvectors.cols()==neig);

  // this function does not tune the density estimation
  const bool tuneDens = false;

  // compute the kernel matrix and the diagonal matrix S^{-1}
  Eigen::SparseMatrix<double> kmat;
  S = KernelMatrix(BandwidthParameter(), kmat, &tuneDens);
  S = P.array()*S.array().sqrt();
  Sinv = S.array().inverse();

  // compute Lhat, which is related to the Kolmogorov operator by a similarity transformation
  Eigen::SparseMatrix<double> Lhat(n, n);
  Lhat = (P.array()*P.array()).inverse().matrix().asDiagonal();
  Lhat = (Sinv.asDiagonal()*kmat*Sinv.asDiagonal()-Lhat)/BandwidthParameter();

  // compute the eigendecomposition of Lhat
  Spectra::SparseGenMatProd<double> matvec(Lhat);
  Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolver(&matvec, neig, std::min(10*neig, n));
  eigsolver.init();
  const std::size_t ncomputed = eigsolver.compute(eigensolverMaxIt, eigensolverTol);
  assert(ncomputed==neig);

  // fill the eigenvectors/values
  eigenvalues = eigsolver.eigenvalues();
  eigenvectors = eigsolver.eigenvectors();
}

std::size_t KolmogorovOperator::NumEigenvalues() const { return neig; }

double KolmogorovOperator::EigensolverTolerance() const { return eigensolverTol; }

double KolmogorovOperator::EigensolverMaxIterations() const { return eigensolverMaxIt; }

Eigen::VectorXd KolmogorovOperator::FunctionRepresentation(Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, std::function<double(Eigen::VectorXd const&)> f) const { return FunctionRepresentation(S.asDiagonal()*eigenvectors, f); }

Eigen::VectorXd KolmogorovOperator::FunctionRepresentation(Eigen::Ref<const Eigen::MatrixXd> const& eigenvectorsRight, std::function<double(Eigen::VectorXd const&)> f) const {
  // the number of samples
  const std::size_t n = NumSamples();

  // evaluate the function at each sample
  Eigen::VectorXd feval(n);
  for( std::size_t i=0; i<n; ++i ) { feval(i) = f(Point(i)); }

  return FunctionRepresentation(eigenvectorsRight, feval);
}

Eigen::VectorXd KolmogorovOperator::FunctionRepresentation(Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, Eigen::Ref<const Eigen::VectorXd> const& feval) const { return FunctionRepresentation(S.asDiagonal()*eigenvectors, feval); }

Eigen::VectorXd KolmogorovOperator::FunctionRepresentation(Eigen::Ref<const Eigen::MatrixXd> const& eigenvectorsRight, Eigen::Ref<const Eigen::VectorXd> const& feval) const {
  assert(feval.size()==NumSamples());
  return eigenvectorsRight.transpose()*feval;
}

Eigen::VectorXd KolmogorovOperator::PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& eigenvalues) const {
  Eigen::VectorXd eigenvaluesInv(eigenvalues.size());
  for( std::size_t i=0; i<eigenvalues.size(); ++i ) {
    eigenvaluesInv(i) = (std::abs(eigenvalues(i))>2.0*eigensolverTol? 1.0/eigenvalues(i) : 0.0);
  }

  return eigenvaluesInv;
}

Eigen::VectorXd KolmogorovOperator::PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& rhs, Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::VectorXd> const& Sinv, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, bool const inv) const {
  if( !inv ) {
    // compute the eigenvalue inverse
    const Eigen::VectorXd eigenvaluesInv = PseudoInverse(eigenvalues);
    return PseudoInverse(rhs, S, Sinv, eigenvaluesInv, eigenvectors, true);
  }

  // compute the pseudo inverse
  Eigen::VectorXd Linv_rhs = Sinv.asDiagonal()*eigenvectors*PseudoInverse(rhs, S, eigenvalues, eigenvectors, true);
  assert(Linv_rhs.size()==NumSamples());

  // make the expected value zero
  Linv_rhs -= Eigen::VectorXd::Constant(Linv_rhs.size(), Linv_rhs.sum()/Linv_rhs.size());

  return Linv_rhs;
}

Eigen::VectorXd KolmogorovOperator::PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& rhs, Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, bool const inv) const {
  assert(rhs.size()==NumSamples());

  if( !inv ) {
    // compute the eigenvalue inverse
    const Eigen::VectorXd eigenvaluesInv = PseudoInverse(eigenvalues);
    return PseudoInverse(rhs, S, eigenvaluesInv, eigenvectors, true);
  }

  // compute the function representation of the rhs
  const Eigen::VectorXd rhsCoeff = FunctionRepresentation(S, eigenvectors, rhs);

  return PseudoInverse(rhsCoeff, eigenvalues, true);
}

Eigen::VectorXd KolmogorovOperator::PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& rhs, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, bool const inv) const {
  assert(rhs.size()==eigenvalues.size());

  if( !inv ) {
    // compute the eigenvalue inverse
    const Eigen::VectorXd eigenvaluesInv = PseudoInverse(eigenvalues);
    return PseudoInverse(rhs, eigenvaluesInv, true);
  }

  return eigenvalues.asDiagonal()*rhs;
}

Eigen::MatrixXd KolmogorovOperator::FunctionGradient(Eigen::Ref<const Eigen::VectorXd> const& coeff, Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::VectorXd> const& Sinv, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors) const {
  // the coordinate dimension and the number of samples
  const std::size_t dim = Point(0).size();
  const std::size_t n = NumSamples();

  // compute the left and right eigenvectors
  const Eigen::MatrixXd eigenvectorsLeft = Sinv.asDiagonal()*eigenvectors;
  const Eigen::MatrixXd eigenvectorsRight = S.asDiagonal()*eigenvectors;

  // precompute the coeffients for each dimension
  Eigen::MatrixXd xcoeff(neig, dim);
  for( std::size_t d=0; d<dim; ++d ) {
    xcoeff.col(d) = FunctionRepresentation(eigenvectorsRight, [d](Eigen::VectorXd const& x) -> double { return x(d); });
  }

  // compute the coordinate coefficients
  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(n, dim);
  for( std::size_t j=0; j<neig; ++j ) {
    for( std::size_t k=j; k<neig; ++k ) {
      const Eigen::VectorXd phijk = eigenvectorsLeft.col(j).array()*eigenvectorsLeft.col(k).array();
      assert(phijk.size()==n);

      for( std::size_t l=0; l<neig; ++l ) {
        const double Cjkl = phijk.dot(eigenvectorsRight.col(l))/2.0;
        const double scale = Cjkl*(eigenvalues(l)-eigenvalues(k)-eigenvalues(j));
        const double scalej = coeff(j)*scale;

        for( std::size_t d=0; d<dim; ++d ) {
          gradient.col(d) += (scalej*xcoeff.col(d)(k))*eigenvectorsLeft.col(l);

          if( j!=k ) {
            const double scalek = coeff(k)*scale;
            for( std::size_t d=0; d<dim; ++d ) {
              gradient.col(d) += (scalek*xcoeff.col(d)(j))*eigenvectorsLeft.col(l);
            }
          }
        }
      }
    }
  }

  return gradient;
}
