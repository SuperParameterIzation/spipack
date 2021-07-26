#include <iostream>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp>

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace spi::NumericalSolvers;

int main(int argc, char **argv) {
  // the dimension of the problem
  const std::size_t dim = 2;

  // the output filename
  const std::string filename = "outputData.h5";

  // the number of samples
  const std::size_t n = 10000;

  // numerical parameters
  const std::size_t numNeighbors = 25;

  // the number of eigenvalues
  const std::size_t neigs = 3*log(n);
  const double eigensolverTol = 1.0e-4;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = std::size_t(5*std::log((double)n));
  nnOptions["NumThreads"] = omp_get_max_threads();

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "ExponentialKernel";

  // set the options for the density estimation
  YAML::Node densityOptions;
  densityOptions["NearestNeighbors"] = nnOptions;
  densityOptions["KernelOptions"] = kernelOptions;
  densityOptions["NumNearestNeighbors"] = numNeighbors;
  densityOptions["TruncationTolerance"] = -std::log(1.0e-4);
  densityOptions["ManifoldDimension"] = (double)dim;

  // set the options for the Kolmogorov operator
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["KernelOptions"] = kernelOptions;
  options["DensityOptions"] = densityOptions;
  options["OperatorParameter"] = 1.0;
  options["TruncationTolerance"] = -std::log(1.0e-4);
  options["ManifoldDimension"] = (double)dim;
  options["BandwidthExponent"] = -0.5;
  options["BandwidthParameter"] = 1.0e-1;
  options["NumEigenvalues"] = neigs;
  options["EigensolverTolerance"] = eigensolverTol;
  options["NumNearestNeighbors"] = 25;

  // create the Kolmogorov operator
  auto kolOperator = std::make_shared<KolmogorovOperator>(rv, options);

  // construct the kd-trees
  std::cout << "building kd trees ... " << std::flush;
  kolOperator->BuildKDTrees();
  std::cout << "done." << std::endl;

  // tune the bandwidth parameter
  std::cout << "tune the bandwidth parameter ... " << std::flush;
  kolOperator->TuneBandwidthParameter(true);
  std::cout << "done." << std::endl;

  // compute the eigendecomposition of Lhat
  std::cout << "computing eigendecomposition of Lhat ... " << std::flush;
  Eigen::VectorXd S(n), Sinv(n), lambda(neigs);
  Eigen::MatrixXd Qhat(n, neigs);
  kolOperator->ComputeEigendecomposition(S, Sinv, lambda, Qhat);
  std::cout << "done." << std::endl;

  // solve the linear system using the eigendecomposition
  std::cout << "computing solution to the inverse problem ... " << std::flush;
  //const Eigen::VectorXd direction = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd direction = Eigen::VectorXd::Zero(dim);
  direction(0) = 1.0;
  Eigen::VectorXd rhs(n);
  for( std::size_t i=0; i<n; ++i ) { rhs(i) = direction.dot(kolOperator->Point(i)); }

  // apply the pseudo-inverse to the rhs (we will need the coeffients of the solution)
  const Eigen::VectorXd inverseCoeff = kolOperator->PseudoInverse(rhs, S, lambda, Qhat);
  assert(inverseCoeff.size()==kolOperator->NumEigenvalues());
  Eigen::VectorXd inverse = Sinv.asDiagonal()*Qhat*inverseCoeff;
  //inverse -= Eigen::VectorXd::Constant(n, inverse.sum()/n);


  const Eigen::VectorXd rhsCoeff = kolOperator->FunctionRepresentation(S, Qhat, rhs);
  Eigen::VectorXd rhsExpanded = Sinv.asDiagonal()*Qhat*rhsCoeff;

  std::cout << "done." << std::endl;

  // solve the linear system using the eigendecomposition
  std::cout << "computing the solution gradient ... " << std::flush;

  const Eigen::MatrixXd gradient = kolOperator->FunctionGradient(inverseCoeff, S, Sinv, lambda, Qhat);

  std::cout << "done." << std::endl;

  // write the samples to file
  kolOperator->WriteToFile(filename);

  // open the file and write data to file
  HDF5File hdf5file(filename);
  hdf5file.WriteMatrix("/true right hand side", rhs);
  hdf5file.WriteMatrix("/expanded right hand side", rhsExpanded);
  hdf5file.WriteMatrix("/weighted Poisson solution", inverse);
  hdf5file.WriteMatrix("/weighted Poisson solution gradient", gradient);
  hdf5file.WriteMatrix("/eigenvectors Qhat", Qhat);
  hdf5file.Close();
}
