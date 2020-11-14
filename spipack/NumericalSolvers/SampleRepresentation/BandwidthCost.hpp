#ifndef BANDWIDTHCOST_HPP_
#define BANDWIDTHCOST_HPP_

#include <MUQ/Optimization/CostFunction.h>

#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

namespace spi {
namespace NumericalSolvers {

class BandwidthCost : public muq::Optimization::CostFunction {
public:

  BandwidthCost(std::shared_ptr<SampleRepresentation> const& samples);

  virtual ~BandwidthCost() = default;
private:
   virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

   virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

   std::shared_ptr<SampleRepresentation> samples;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
