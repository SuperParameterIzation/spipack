#include <iostream>

#include <MUQ/Modeling/Distributions/UniformBox.h>

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

  // create a uniform random variable
  std::vector<std::pair<double, double> > bounds(dim, std::pair<double, double>(1.0, 0.0));
  auto rv = std::make_shared<UniformBox>(bounds)->AsVariable();

  // the options for the graph Laplacian
  YAML::Node options;
  options["NumSamples"] = 1000;
  options["Bandwidth"] = 0.25;
  options["EigensolverTol"] = 1.0e-4;

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "HatKernel";
  options["KernelOptions"] = kernelOptions;

  // create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(rv, options);

  // write the samples to file
  laplacian->WriteToFile(filename);

  // construct the heat matrix
  laplacian->ConstructHeatMatrix();

  // the largest eigenvalue is 1
  const std::size_t neig = 100;
  const Eigen::VectorXd eigenvalues = laplacian->HeatMatrixEigenvalues(neig);

  // open the file
  auto hdf5file = std::make_shared<HDF5File>(filename);
  hdf5file->WriteMatrix("/heat matrix eigenvalues", eigenvalues);
  hdf5file->Close();
}
