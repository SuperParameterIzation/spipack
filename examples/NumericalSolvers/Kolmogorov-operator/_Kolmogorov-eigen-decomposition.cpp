#include <iostream>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <MUQ/Modeling/Distributions/Density.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/Tools/SortVector.hpp>

#include <spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp>
#include <spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp>

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;

int main(int argc, char **argv) {
  // the dimension of the problem
  const std::size_t dim = 2;

  // the output filename
  const std::string filename = "outputData_eigendecomposition.h5";

  // the number of samples
  const std::size_t n = 10000;

  // numerical parameters
  const std::size_t numNeighbors = 25;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // compute a sample collection
  auto samples = std::make_shared<SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { samples->Add(std::make_shared<SamplingState>(rv->Sample())); }

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = omp_get_max_threads();

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "ExponentialKernel";

  // set the options for the density estimation
  YAML::Node densityOptions;
  densityOptions["NearestNeighbors"] = nnOptions;
  densityOptions["KernelOptions"] = kernelOptions;
  densityOptions["NumNearestNeighbors"] = numNeighbors;
  densityOptions["TruncationTolerance"] = -std::log(1.0e-2);
  densityOptions["ManifoldDimension"] = (double)dim;

  // set the options for the Kolmogorov operator
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["KernelOptions"] = kernelOptions;
  options["DensityOptions"] = densityOptions;
  options["OperatorParameter"] = 1.0;
  options["TruncationTolerance"] = -std::log(1.0e-2);
  options["ManifoldDimension"] = (double)dim;
  options["BandwidthExponent"] = -0.5;
  options["BandwidthParameter"] = 1.0e-1;

  // create the Kolmogorov operator
  auto kolOperator = std::make_shared<KolmogorovOperator>(samples, options);

  // construct the kd-trees
  std::cout << "building kd trees ... " << std::flush;
  kolOperator->BuildKDTrees();
  std::cout << "done." << std::endl;

  // tune the bandwidth parameter
  std::cout << "tune the bandwidth parameter ... " << std::flush;
  kolOperator->TuneBandwidthParameter(true);
  std::cout << "done." << std::endl;

  // estimate the density at each sample
  std::cout << "estimating the underlying density ... " << std::flush;
  Eigen::VectorXd dens = kolOperator->EstimateDensity(false);
  std::cout << "done." << std::endl;

  // compute the tuning cost function (for plotting purposes)
  std::cout << "computing tuning cost function (for plotting purposes only) ... " << std::flush;
  BandwidthCost bandCost(dens.array().pow(kolOperator->ExponentParameter()), kolOperator->Density());
  const Eigen::VectorXd bandwidthExponent = Eigen::VectorXd::LinSpaced(50, -25.0, 5.0);
  Eigen::VectorXd candidateBandwidthParameters(bandwidthExponent.size());
  Eigen::VectorXd logKernelAvgChange(bandwidthExponent.size());
  for( std::size_t l=0; l<bandwidthExponent.size(); ++l ) {
    candidateBandwidthParameters(l) = std::pow(2.0, bandwidthExponent(l))/4.0;
    logKernelAvgChange(l) = -bandCost.Cost(Eigen::VectorXd::Constant(1, bandwidthExponent(l)).eval());
  }
  std::cout << "done." << std::endl;

  // estimate the symmetric kernel matrix K
  std::cout << "estimating the kernel matrix K ... " << std::flush;
  Eigen::SparseMatrix<double> kmat;
  const Eigen::VectorXd D = kolOperator->KernelMatrix(kolOperator->BandwidthParameter(), dens, kmat);
  std::cout << "done." << std::endl;

  // compute the diagonal matrices P=\psi^{beta} and S = P D^{1/2}
  std::cout << "creating L and Lhat ... " << std::flush;
  const Eigen::VectorXd P = dens.array().pow(kolOperator->ExponentParameter());
  const Eigen::VectorXd S = P.array()*(D.array().sqrt());
  const Eigen::VectorXd Sinv = S.array().inverse();

  // directly compute the discrete Kolmogorov operator
  Eigen::SparseMatrix<double> Lkol(n, n);
  Lkol.setIdentity();
  Lkol = (P.array()*P.array()).inverse().matrix().asDiagonal()*(D.array().inverse().matrix().asDiagonal()*kmat-Lkol)/kolOperator->BandwidthParameter();

  // a test function
  Eigen::MatrixXd f = Eigen::MatrixXd::Ones(n, 3); // f(x) = 1
  for( std::size_t i=0; i<n; ++i ) {
    f(i, 1) = samples->at(i)->state[0](0); // f(x) = x_0
    f(i, 2) = samples->at(i)->state[0](1); // f(x) = x_1
  }

  // apply the discrete Kolmogorov operator to a function
  const Eigen::MatrixXd Lf = Lkol*f;

  // compute Lhat, which is related to the Kolmogorov operator by a similarity transformation
  Eigen::SparseMatrix<double> Lhat(n, n);
  Lhat = (P.array()*P.array()).inverse().matrix().asDiagonal();
  Lhat = (Sinv.asDiagonal()*kmat*Sinv.asDiagonal()-Lhat)/kolOperator->BandwidthParameter();

  // apply the similarity trasformation operator to a function
  const Eigen::MatrixXd SinvLhatSf = Sinv.asDiagonal()*Lhat*S.asDiagonal()*f;
  std::cout << "done." << std::endl;

  // compute the eigendecomposition of Lhat
  std::cout << "computing eigendecomposition of Lhat ... " << std::flush;
  Spectra::SparseGenMatProd<double> matvec(Lhat);
  const std::size_t neigs = 100;
  Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolver(&matvec, neigs, std::min(10*neigs, n));
  eigsolver.init();
  const std::size_t ncomputed = eigsolver.compute(1000, 1.0e-5);
  std::cout << "done (" << ncomputed << " of " << neigs << " computed)" << std::endl;

  // retrieve the eigenvalues and eigen vectors (of Lhat)
  const Eigen::VectorXd lambda = eigsolver.eigenvalues();
  const Eigen::MatrixXd Uhat = eigsolver.eigenvectors();
  Eigen::VectorXd lambdaInv(neigs);
  for( std::size_t i=0; i<neigs; ++i ) { lambdaInv(i) = (std::abs(lambda(i))>1.0e-12? 1.0/lambda(i) : 0.0); }

  // apply Lhat inverse to the function
  const Eigen::MatrixXd Lhatinvf = Uhat*lambdaInv.asDiagonal()*Uhat.transpose()*f;

  // apply L inverse to the function
  const Eigen::MatrixXd Linvf = Sinv.asDiagonal()*Uhat*lambdaInv.asDiagonal()*Uhat.transpose()*S.asDiagonal()*f;

  // write the samples to file
  kolOperator->WriteToFile(filename);

    // open the file and write data to file
  HDF5File hdf5file(filename);
  hdf5file.WriteMatrix("/density estimate", dens);
  hdf5file.WriteMatrix("/eigenvalues", lambda);
  hdf5file.WriteMatrix("/eigenvectors (L)", (Sinv.asDiagonal()*Uhat).eval());
  hdf5file.WriteMatrix("/eigenvectors (Lhat)", Uhat);

  hdf5file.WriteMatrix("/tune/candidate bandwidth parameters", candidateBandwidthParameters);
  hdf5file.WriteMatrix("/tune/log kernel average change", logKernelAvgChange);

  hdf5file.WriteMatrix("/apply function", f);
  hdf5file.WriteMatrix("/applied Kolmogorov operator (L)", Lf);
  hdf5file.WriteMatrix("/applied Kolmogorov operator (Sinv Lhat S)", SinvLhatSf);
  hdf5file.WriteMatrix("/applied inverse transfromed Kolmogorov operator", Lhatinvf);
  hdf5file.WriteMatrix("/applied inverse Kolmogorov operator", Linvf);
  hdf5file.Close();
}
