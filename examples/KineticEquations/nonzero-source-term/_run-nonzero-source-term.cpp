#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/KineticEquations/ConditionalVelocityDistribution.hpp>

using namespace muq::Modeling;
using namespace spi::KineticEquations;

// the dimension of the state space
constexpr std::size_t dim = 2;

std::shared_ptr<ConditionalVelocityDistribution::MacroscaleInformation> MacroscaleInformation() {
  Eigen::VectorXd logMassDensityGrad = Eigen::Matrix<double, dim, 1>::Zero(dim, 1);
  logMassDensityGrad(0) = 1.0;

  return std::make_shared<ConditionalVelocityDistribution::MacroscaleInformation>(1.0, Eigen::Matrix<double, dim, 1>::Zero(dim, 1), 0.0, logMassDensityGrad);
  }

int main(int argc, char **argv) {
  // the number of samples
  const std::size_t n = 1000;

  // the number of timesteps
  const std::size_t numTimesteps = 100;

  // the output filename
  const std::string filename = "output/output";

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = omp_get_max_threads();

  // options for the parameter tuning optimization
  YAML::Node optOptions;
  optOptions["Algorithm"] = "SBPLX";
  optOptions["Ftol.AbsoluteTolerance"] = 1.0e-1;
  optOptions["Ftol.RelativeTolerance"] = 1.0e-1;
  optOptions["Rtol.AbsoluteTolerance"] = 1.0e-1;
  optOptions["Rtol.RelativeTolerance"] = 1.0e-1;

  // options for the kolmogorov operator
  YAML::Node kolOptions;
  kolOptions["Optimization"] = optOptions;
  kolOptions["NumNearestNeighbors"] = 50;
  kolOptions["TruncationTolerance"] = -std::log(1.0e-4);

  // set the options for the conditional velocity distribution
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["NumTimesteps"] = numTimesteps;
  options["KolmogorovOptions"] = kolOptions;
  options["OutputFilename"] = filename;
  options["AccelerationNoiseScale"] = 0.0;
  options["KolmogorovRetuneFrequency"] = 1000;

  auto distribution = std::make_shared<ConditionalVelocityDistribution>(
    Eigen::Matrix<double, dim, 1>::Zero(dim),
    std::make_shared<Gaussian>(dim)->AsVariable(), MacroscaleInformation(),
    options);

  // run the micro-scale model
  double time = 0.0;
  const double delta = 0.05;
  for( std::size_t i=0; i<5; ++i ) {
    time += delta;
    distribution->Run(time, MacroscaleInformation(), i==0);
  }
}
