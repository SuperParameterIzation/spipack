#include <iostream>

#include <MUQ/Modeling/Distributions/Density.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp>

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
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

  // create a density/random variable
  auto gauss = std::make_shared<Gaussian>(dim);
  auto logDensity = gauss->AsDensity();
  auto rv = gauss->AsVariable();

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

  // the options for the graph Laplacian
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["KernelOptions"] = kernelOptions;
  options["NumNearestNeighbors"] = numNeighbors;
  options["TruncationTolerance"] = -std::log(1.0e-1);
  options["ManifoldDimension"] = (double)dim;

  // create the density estimator
  DensityEstimation density(samples, options);

  // build the kd trees
  density.BuildKDTrees();

  // compute the squared bandwidth parameter
  const Eigen::VectorXd squaredBandwidth = density.SquaredBandwidth();

  // create the tuning data
  DensityEstimation::TuningData tune;
  tune.bandwidthExponent = Eigen::VectorXd::LinSpaced(50, -10.0, 10.0);

  // estimate the density at each sample
  const Eigen::VectorXd dens = density.Estimate(squaredBandwidth, tune);

  // estimate the density at each sample
  //const Eigen::VectorXd dens = density.Estimate(squaredBandwidth);

  // compute the true density
  Eigen::VectorXd logDens(n);
  for( std::size_t i=0; i<n; ++i ) { logDens(i) = logDensity->LogDensity(samples->at(i)->state[0]); }

  // write the samples to file
  density.WriteToFile(filename);

  // open the file and write data to file
  HDF5File hdf5file(filename);
  hdf5file.WriteMatrix("/squared bandwidth", squaredBandwidth);
  hdf5file.WriteMatrix("/density estimate", dens);
  hdf5file.WriteMatrix("/true density", logDens.array().exp().matrix().eval());

  hdf5file.WriteMatrix("/tune/candidate bandwidth parameters", tune.candidateBandwidthParameters);
  hdf5file.WriteMatrix("/tune/kernel average", tune.kernelAvg);
  hdf5file.WriteMatrix("/tune/log kernel average derivative", tune.logKernelAvgDerivative);
  hdf5file.Close();

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
