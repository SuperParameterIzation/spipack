#include "spipack/NumericalSolvers/SampleRepresentation/BandwidthCost.hpp"

using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace spi::NumericalSolvers;

BandwidthCost::BandwidthCost(Eigen::Ref<const Eigen::VectorXd> const& bandwidth, std::shared_ptr<SampleRepresentation> const& samples) :
CostFunction(Eigen::VectorXi::Constant(1, 1)),
bandwidth(bandwidth),
samples(samples)
{}

double BandwidthCost::CostImpl(ref_vector<Eigen::VectorXd> const& input) {
  const double neweps = std::pow(2.0, input[0](0));
  const std::size_t n2 = samples->NumSamples()*samples->NumSamples();

  if( std::abs(eps-neweps)>1.0e-10 || std::isnan(eps) ) {
    eps = neweps;
    Eigen::SparseMatrix<double> kmat;
    Eigen::VectorXd rowsum = samples->KernelMatrix(std::pow(2.0, input[0](0)+1.0), bandwidth, kmat);
    ksumlp1 = rowsum.sum()/n2;

    rowsum = samples->KernelMatrix(eps, bandwidth, kmat);

    ksuml = rowsum.sum()/n2;
  }

  const double siglp1 = std::log(ksumlp1);
  const double sigl = std::log(ksuml);
  const double cost = (siglp1-sigl)/std::log(2.0);

  return -cost;
}

void BandwidthCost::GradientImpl(unsigned int const inputDimWrt, ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  assert(inputDimWrt==0);
  const double neweps = std::pow(2.0, input[0](0));
  const double epsp1 = 2.0*neweps;
  const std::size_t n2 = samples->NumSamples()*samples->NumSamples();

  if( std::abs(eps-neweps)>1.0e-10 || std::isnan(eps) ) {
    eps = neweps;

    Eigen::SparseMatrix<double> kmat;
    Eigen::VectorXd rowsum = samples->KernelMatrix(epsp1, bandwidth, kmat);

    ksumlp1 = rowsum.sum()/n2;

    rowsum = samples->KernelMatrix(eps, bandwidth, kmat);

    ksuml = rowsum.sum()/n2;
  }

  const double avgl = samples->KernelDerivativeAverage(eps, bandwidth);
  const double avglp1 = samples->KernelDerivativeAverage(epsp1, bandwidth);

  gradient.resize(1);
  gradient(0) = -(avglp1/ksumlp1*epsp1 - avgl/ksuml*eps)/std::log(2.0);
}
