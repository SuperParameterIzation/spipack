#ifndef BANDWIDTHCOST_HPP_
#define BANDWIDTHCOST_HPP_

#include <MUQ/Optimization/CostFunction.h>

namespace spi {
namespace NumericalSolvers {

class BandwidthCost {
public:

  BandwidthCost();

  virtual ~BandwidthCost() = default;
private:
};

} // namespace NumericalSolvers
} // namespace spi

#endif
