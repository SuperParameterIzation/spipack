#include <iostream>

#include <Spectra/SymEigsSolver.h>
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
  const std::string filename = "outputData.h5";

  // the number of samples
  const std::size_t n = 1000;

  // numerical parameters
  const std::size_t numNeighbors = 25;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

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
  densityOptions["TruncationTolerance"] = -std::log(1.0e-1);
  densityOptions["ManifoldDimension"] = (double)dim;

  // set the options for the Kolmogorov operator
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["KernelOptions"] = kernelOptions;
  options["DensityOptions"] = densityOptions;
  options["OperatorParameter"] = 1.0;
  options["TruncationTolerance"] = -std::log(1.0e-1);
  options["ManifoldDimension"] = (double)dim;

  // create the Kolmogorov operator
  auto kolOperator = std::make_shared<KolmogorovOperator>(rv, options);

  // construct the kd-trees
  kolOperator->BuildKDTrees();

  // tune the bandwidth parameter
  kolOperator->TuneBandwidthParameter();

  // estimate the density at each sample
  Eigen::VectorXd dens = kolOperator->EstimateDensity(false);

  //std::cout << eps*(dens.array()*dens.array()).matrix().transpose() << std::endl;

  Eigen::SparseMatrix<double> kmat;
  Eigen::VectorXd D = kolOperator->KernelMatrix(kolOperator->BandwidthParameter(), dens, kmat);

  const Eigen::VectorXd P = dens.array().pow(kolOperator->ExponentParameter());
  const Eigen::VectorXd S = P.array()*(D.array().sqrt());

  // create a cost function
  BandwidthCost bandCost(P, kolOperator);

  // create the tuning data
  const Eigen::VectorXd bandwidthExponent = Eigen::VectorXd::LinSpaced(50, -10.0, 10.0);
  Eigen::VectorXd candidateBandwidthParameters(bandwidthExponent.size());
  Eigen::VectorXd logKernelAvgChange(bandwidthExponent.size());
  for( std::size_t l=0; l<bandwidthExponent.size(); ++l ) {
    candidateBandwidthParameters(l) = std::pow(2.0, bandwidthExponent(l));
    logKernelAvgChange(l) = -bandCost.Cost(Eigen::VectorXd::Constant(1, bandwidthExponent(l)).eval());
  }

  //Eigen::MatrixXd Lhat = (dens.array()*dens.array()).inverse().matrix().asDiagonal();
  Eigen::SparseMatrix<double> Lhat(n, n);
  Lhat.setIdentity();

  //Lhat = (P.array()*P.array()).inverse().matrix().asDiagonal();

  //Lhat = (S.array().inverse().matrix().asDiagonal()*kmat*S.array().inverse().matrix().asDiagonal() - Lhat)/kolOperator->BandwidthParameter();
  Lhat = (D.array().sqrt().inverse().matrix().asDiagonal()*kmat*D.array().sqrt().inverse().matrix().asDiagonal() - Lhat)/kolOperator->BandwidthParameter();
  //Lhat = D.asDiagonal()*kmat - Lhat;///(eps*eps);

  Spectra::SparseGenMatProd<double> matvec(Lhat);

  const std::size_t neigs = 25;

  std::cout << "computing eigs... " << std::endl;

  Spectra::SymEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolver(&matvec, neigs, 10*neigs);

  // initialize and compute
  eigsolver.init();
  const int ncomputed = eigsolver.compute(1000, 1.0e-4);

  std::cout << "ncomputed: " << ncomputed << std::endl;

  Eigen::VectorXd eigvals = eigsolver.eigenvalues();
  //Eigen::MatrixXd eigvecs = (dens.array().pow(-2.0*kolOperator->ExponentParameter())*D.array().sqrt().inverse()).matrix().asDiagonal()*eigsolver.eigenvectors()/kolOperator->BandwidthParameter();
  //Eigen::MatrixXd eigvecs = eigsolver.eigenvectors();
  Eigen::MatrixXd eigvecs = eigsolver.eigenvectors();


  //std::cout << eigvecs.transpose()*eigvecs << std::endl;

  Eigen::VectorXd R = Eigen::VectorXd::Ones(n);


  //std::cout << dens.array().pow(2.0).matrix().asDiagonal()*R/eps << std::endl;

  // write the samples to file
  kolOperator->WriteToFile(filename);

    // open the file and write data to file
  HDF5File hdf5file(filename);
  hdf5file.WriteMatrix("/density estimate", dens);
  hdf5file.WriteMatrix("/eigenvalues", eigvals);
  hdf5file.WriteMatrix("/eigenvectors", eigvecs);

  hdf5file.WriteMatrix("/tune/candidate bandwidth parameters", candidateBandwidthParameters);
  hdf5file.WriteMatrix("/tune/log kernel average change", logKernelAvgChange);
  hdf5file.Close();
  hdf5file.Close();
}
