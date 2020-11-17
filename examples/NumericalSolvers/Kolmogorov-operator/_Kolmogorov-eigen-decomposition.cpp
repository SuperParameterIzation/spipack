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
  const std::string filename = "outputData.h5";

  // the number of samples
  const std::size_t n = 1000;

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
  auto kolOperator = std::make_shared<KolmogorovOperator>(samples, options);

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

  std::cout << "finished set up" << std::endl;

  // create the discrete Kolmogorov operator
  Eigen::SparseMatrix<double> Lkol(n, n);
  Lkol.setIdentity();
  Lkol = (P.array()*P.array()).inverse().matrix().asDiagonal()*(D.array().inverse().matrix().asDiagonal()*kmat-Lkol)/kolOperator->BandwidthParameter();

  // application to a function
  Eigen::MatrixXd f = Eigen::MatrixXd::Ones(n, 3);
  for( std::size_t i=0; i<n; ++i ) {
    f(i, 1) = samples->at(i)->state[0](0);
    f(i, 2) = samples->at(i)->state[0](1);
  }
  const Eigen::MatrixXd Lf = Lkol*f;

  /*// create a copy of the discrete operator that enforces E[h] = 0
  Eigen::SparseMatrix<double> Lboundary = Lkol;
  for( std::size_t i=0; i<n; ++i ) {
    Lboundary.coeffRef(0, i) = 1.0;
  }

  // invert the discrete laplacian (enforce that the expected value is 0)
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
  solver.compute(Lboundary);
  Eigen::MatrixXd rhs = f;
  rhs.row(0) = Eigen::VectorXd::Zero(3);
  Eigen::MatrixXd Linvf = solver.solve(rhs);*/

  //Eigen::MatrixXd Lhat = (dens.array()*dens.array()).inverse().matrix().asDiagonal();
  Eigen::SparseMatrix<double> Lhat(n, n);
  Lhat.setIdentity();
  Lhat = (P.array()*P.array()).inverse().matrix().asDiagonal();
  Lhat = S.array().inverse().matrix().asDiagonal()*kmat*S.array().inverse().matrix().asDiagonal()-Lhat;

  //Lhat = (P.array()*P.array()).inverse().matrix().asDiagonal();

  //Lhat = (S.array().inverse().matrix().asDiagonal()*kmat*S.array().inverse().matrix().asDiagonal() - Lhat)/kolOperator->BandwidthParameter();
  //Lhat = (D.array().sqrt().inverse().matrix().asDiagonal()*kmat*D.array().sqrt().inverse().matrix().asDiagonal() - Lhat)/kolOperator->BandwidthParameter();
  //Lhat = D.asDiagonal()*kmat - Lhat;///(eps*eps);

  Spectra::SparseGenMatProd<double> matvec(Lhat);

  const std::size_t neigs = 25;

  std::cout << "computing eigs... " << std::endl;

  Spectra::GenEigsSolver<double, Spectra::SMALLEST_MAGN, Spectra::SparseGenMatProd<double> > eigsolver(&matvec, neigs, 10*neigs);

  // initialize and compute
  eigsolver.init();
  const int ncomputed = eigsolver.compute(1000, 1.0e-4);

  std::cout << "ncomputed: " << ncomputed << std::endl;

  Eigen::VectorXd eigvals = eigsolver.eigenvalues().real();
  //Eigen::MatrixXd eigvecs = (dens.array().pow(-2.0*kolOperator->ExponentParameter())*D.array().sqrt().inverse()).matrix().asDiagonal()*eigsolver.eigenvectors()/kolOperator->BandwidthParameter();
  //Eigen::MatrixXd eigvecs = eigsolver.eigenvectors();
  Eigen::MatrixXd eigvecs = eigsolver.eigenvectors().real();

  Eigen::VectorXd eigvalsInv(neigs);
  for( std::size_t i=0; i<neigs; ++i ) { eigvalsInv(i) = (std::abs(eigvals(i))>1.0e-12? 1.0/eigvals(i) : 0.0); }

  //Eigen::MatrixXd id = eigvecs*eigvalsInv.asDiagonal()*eigvecs.transpose()*Lhat;

  Eigen::MatrixXd Linvf = kolOperator->BandwidthParameter()*eigvecs*eigvalsInv.asDiagonal()*eigvecs.transpose()*(P.array()*P.array()).matrix().asDiagonal()*f;

  Eigen::MatrixXd Lfeig = (1.0/kolOperator->BandwidthParameter())*(P.array()*P.array()).inverse().matrix().asDiagonal()*eigvecs*eigvals.asDiagonal()*eigvecs.transpose()*f;

  std::cout << eigvecs*eigvals.asDiagonal()*eigvecs.transpose()-Lhat << std::endl;

  //std::cout << id << std::endl;


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

  hdf5file.WriteMatrix("/apply function", f);
  hdf5file.WriteMatrix("/applied Kolmogorov operator", Lf);
  hdf5file.WriteMatrix("/applied Kolmogorov operator (eigendecomposition)", Lfeig);
  hdf5file.WriteMatrix("/applied inverse Kolmogorov operator", Linvf);
  hdf5file.Close();
}
