#include <iostream>

#include <MUQ/Modeling/Distributions/Density.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp>
#include <spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp>

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
  const std::size_t numNeighbors = 50;

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
  auto density = std::make_shared<DensityEstimation>(samples, options);

  // build the kd trees
  density->BuildKDTrees();

  // compute the squared bandwidth parameter
  const Eigen::VectorXd squaredBandwidth = density->SquaredBandwidth();

  // estimate the density at each sample
  const Eigen::VectorXd dens = density->Estimate(squaredBandwidth);

  // create a cost function
  BandwidthCost bandCost(squaredBandwidth.array().sqrt(), density);

  // create the tuning data
  const Eigen::VectorXd bandwidthExponent = Eigen::VectorXd::LinSpaced(50, -10.0, 10.0);
  Eigen::VectorXd candidateBandwidthParameters(bandwidthExponent.size());
  Eigen::VectorXd logKernelAvgChange(bandwidthExponent.size());
  for( std::size_t l=0; l<bandwidthExponent.size(); ++l ) {
    candidateBandwidthParameters(l) = std::pow(2.0, bandwidthExponent(l));
    logKernelAvgChange(l) = -bandCost.Cost(Eigen::VectorXd::Constant(1, bandwidthExponent(l)).eval());
  }

  // compute the true density
  Eigen::VectorXd logDens(n);
  for( std::size_t i=0; i<n; ++i ) { logDens(i) = logDensity->LogDensity(samples->at(i)->state[0]); }

  // write the samples to file
  density->WriteToFile(filename);

  // open the file and write data to file
  HDF5File hdf5file(filename);
  hdf5file.WriteMatrix("/squared bandwidth", squaredBandwidth);
  hdf5file.WriteMatrix("/density estimate", dens);
  hdf5file.WriteMatrix("/true density", logDens.array().exp().matrix().eval());

  hdf5file.WriteMatrix("/tune/candidate bandwidth parameters", candidateBandwidthParameters);
  hdf5file.WriteMatrix("/tune/log kernel average change", logKernelAvgChange);
  hdf5file.Close();
}
