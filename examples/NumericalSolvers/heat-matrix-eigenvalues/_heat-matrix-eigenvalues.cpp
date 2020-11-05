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

  // numerical parameters
  const std::size_t numNeighbors = 10;

  // create a uniform random variable
  std::vector<std::pair<double, double> > bounds(dim, std::pair<double, double>(1.0, 0.0));
  auto rv = std::make_shared<UniformBox>(bounds)->AsVariable();

  // the options for the graph Laplacian
  YAML::Node options;
  options["NumSamples"] = 5000;
  options["Bandwidth"] = 0.25;
  options["EigensolverTol"] = 1.0e-10;

  // set the kernel options
  YAML::Node kernelOptions;
  kernelOptions["Kernel"] = "HatKernel";
  options["KernelOptions"] = kernelOptions;

  // create the graph laplacian
  auto laplacian = std::make_shared<GraphLaplacian>(rv, options);

  // write the samples to file
  laplacian->WriteToFile(filename);

  // build the kd tree
  laplacian->BuildKDTree();

  // loop through each sample
  std::vector<std::vector<std::pair<std::size_t, double> > > neighbors(laplacian->NumSamples());
  for( unsigned int i=0; i<laplacian->NumSamples(); ++i ) {
    // get a reference to the ith point
    Eigen::Ref<Eigen::VectorXd const> point = laplacian->Point(i);

    // find the nearest neighbors for each sample
    //const double h2 = laplacian->FindNeighbors(neighbors[i]);
  }

  /*// construct the heat matrix
  laplacian->ConstructHeatMatrix();

  // the largest eigenvalue is 1
  const std::size_t neig = 100;
  const Eigen::VectorXd eigenvalues = laplacian->HeatMatrixEigenvalues(neig);

  // open the file
  auto hdf5file = std::make_shared<HDF5File>(filename);
  hdf5file->WriteMatrix("/heat matrix eigenvalues", eigenvalues);
  hdf5file->Close();*/
}
