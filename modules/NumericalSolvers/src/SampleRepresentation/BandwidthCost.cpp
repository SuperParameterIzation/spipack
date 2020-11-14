#include "spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace spi::NumericalSolvers;

BandwidthCost::BandwidthCost(std::shared_ptr<SampleRepresentation> const& samples) :
CostFunction(Eigen::Vector2i(1, samples->NumSamples())),
samples(samples)
{}

double BandwidthCost::CostImpl(ref_vector<Eigen::VectorXd> const& input) {
  // compute the average kernel derivative
  return samples->KernelDerivativeAverage(input[0](0), input[1]);
}

void BandwidthCost::GradientImpl(unsigned int const inputDimWrt, ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  assert(inputDimWrt==0);

  gradient.resize(1);
  gradient(0) = samples->KernelSecondDerivativeAverage(input[0](0), input[1]);
}
