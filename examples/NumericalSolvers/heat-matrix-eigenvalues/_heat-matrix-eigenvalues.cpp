#include <iostream>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/NumericalSolvers/GraphLaplacian/GraphLaplacian.hpp>

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::NumericalSolvers;

int main(int argc, char **argv) {
  // the dimension of the problem
  const std::size_t dim = 2;

  // the output filename
  const std::string filename = "samples.h5";

  // the number of samples
  const std::size_t n = 10000;

  // numerical parameters
  const std::size_t numNeighbors = 25.0*std::log(n);
  //const std::size_t numNeighbors = 10;

  // create a uniform random variable
  //std::vector<std::pair<double, double> > bounds(dim, std::pair<double, double>(1.0, 0.0));
  //auto rv = std::make_shared<UniformBox>(bounds)->AsVariable();
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = omp_get_max_threads();

  // the options for the graph Laplacian
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["NumNearestNeighbors"] = numNeighbors;
  options["BandwidthParameter"] = 1.0;
  options["ManifoldDimension"] = 2.0;

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "BumpKernel";
  options["KernelOptions"] = kernelOptions;

  /*// create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(rv, options);

  // write the samples to file
  laplacian->WriteToFile(filename);

  // build the kd tree
  laplacian->BuildKDTrees();

  // compute the bandwidth
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors;
  const Eigen::VectorXd squaredBandwidth = laplacian->SquaredBandwidth(neighbors);

  // compute sigma prime
  std::pair<double, Eigen::SparseMatrix<double> > optimalApprox;
  const Eigen::Matrix<double, Eigen::Dynamic, 2> sigmaprimeApprox = laplacian->TuneKernelBandwidth(squaredBandwidth.array().sqrt(), neighbors, optimalApprox);

  const Eigen::VectorXd densityEstimationApprox = (1.0/(optimalApprox.first*M_PI*n))*optimalApprox.second*squaredBandwidth.array().inverse().matrix();

  std::pair<double, Eigen::SparseMatrix<double> > optimal;
  const Eigen::Matrix<double, Eigen::Dynamic, 2> sigmaprime = laplacian->TuneKernelBandwidth(squaredBandwidth.array().sqrt(), optimal);

  const Eigen::VectorXd densityEstimation = (1.0/(optimal.first*M_PI*n))*optimal.second*squaredBandwidth.array().inverse().matrix();

  // open the file
  auto hdf5file = std::make_shared<HDF5File>(filename);
  //hdf5file->WriteMatrix("/heat matrix eigenvalues", eigenvalues);
  hdf5file->WriteMatrix("/log bandwidth", squaredBandwidth.array().sqrt().log().matrix().eval());
  hdf5file->WriteMatrix("/exact/bandwidth parameter candidate", sigmaprime.col(0).eval());
  hdf5file->WriteMatrix("/exact/sigma prime", sigmaprime.col(1).eval());
  hdf5file->WriteMatrix("/exact/density estimation", densityEstimation);
  hdf5file->WriteMatrix("/approx/bandwidth parameter candidate", sigmaprimeApprox.col(0).eval());
  hdf5file->WriteMatrix("/approx/sigma prime", sigmaprimeApprox.col(1).eval());
  hdf5file->WriteMatrix("/approx/density estimation", densityEstimationApprox);
  hdf5file->Close();*/
}
