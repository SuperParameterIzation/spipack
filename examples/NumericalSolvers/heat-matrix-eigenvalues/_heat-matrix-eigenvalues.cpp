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

  // numerical parameters
  const std::size_t numNeighbors = 10;

  // create a uniform random variable
  //std::vector<std::pair<double, double> > bounds(dim, std::pair<double, double>(1.0, 0.0));
  //auto rv = std::make_shared<UniformBox>(bounds)->AsVariable();
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // the options for the graph Laplacian
  YAML::Node options;
  options["NumSamples"] = 5000;
  options["BandwidthIndex"] = 4;
  options["EigensolverTol"] = 1.0e-10;
  options["BandwidthRange.Min"] = -20.0;
  options["BandwidthRange.Max"] = 10.0;
  options["NumBandwidthSteps"] = 100;

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "BumpKernel";
  options["KernelOptions"] = kernelOptions;

  // create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(rv, options);

  // write the samples to file
  laplacian->WriteToFile(filename);

  // build the kd tree
  laplacian->BuildKDTree();

  // loop through each sample
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors(laplacian->NumSamples());
  Eigen::VectorXd squaredBandwidth(laplacian->NumSamples());
  for( std::size_t i=0; i<laplacian->NumSamples(); ++i ) {
    // get a reference to the ith point
    Eigen::Ref<Eigen::VectorXd const> point = laplacian->Point(i);

    // find the nearest neighbors for each sample
    squaredBandwidth(i) = laplacian->FindNeighbors(point, numNeighbors, neighbors[i]);
  }

  const Eigen::Matrix<double, Eigen::Dynamic, 2> sigmaprime = laplacian->EvaluateKernel(squaredBandwidth.array().sqrt());

  /*for( std::size_t i=0; i<laplacian->NumSamples(); ++i ) {
    // get a reference to the ith point
    Eigen::Ref<Eigen::VectorXd const> point = laplacian->Point(i);

    // evaluate the kernel at this point
    Eigen::MatrixXd kernelEval(numNeighbors+1, laplacian->NumBandwidthSteps()+1);
    laplacian->EvaluateKernel(i, point, neighbors[i], squaredBandwidth.array().sqrt(), kernelEval);
  }*/

  /*// construct the heat matrix
  laplacian->ConstructHeatMatrix();

  // the largest eigenvalue is 1
  const std::size_t neig = 100;
  const Eigen::VectorXd eigenvalues = laplacian->HeatMatrixEigenvalues(neig);*/

  // open the file
  auto hdf5file = std::make_shared<HDF5File>(filename);
  //hdf5file->WriteMatrix("/heat matrix eigenvalues", eigenvalues);
  hdf5file->WriteMatrix("/log squared bandwidth", squaredBandwidth.array().log().matrix().eval());
  hdf5file->WriteMatrix("/bandwidth parameter candidate", sigmaprime.col(0).eval());
  hdf5file->WriteMatrix("/sigma prime", sigmaprime.col(1).eval());
  hdf5file->Close();
}
