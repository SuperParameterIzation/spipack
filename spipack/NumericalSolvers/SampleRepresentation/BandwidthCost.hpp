#ifndef BANDWIDTHCOST_HPP_
#define BANDWIDTHCOST_HPP_

#include <MUQ/Optimization/CostFunction.h>

#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

namespace spi {
namespace NumericalSolvers {

class BandwidthCost : public muq::Optimization::CostFunction {
public:

  BandwidthCost(Eigen::Ref<const Eigen::VectorXd> const& bandwidth, std::shared_ptr<const SampleRepresentation> const& samples);

  virtual ~BandwidthCost() = default;
private:
   virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

   virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

   const Eigen::VectorXd bandwidth;

   std::shared_ptr<const SampleRepresentation> samples;

   double eps = std::numeric_limits<double>::quiet_NaN();
   double ksuml = std::numeric_limits<double>::quiet_NaN();
   double ksumlp1 = std::numeric_limits<double>::quiet_NaN();
};

} // namespace NumericalSolvers
} // namespace spi

#endif
