#ifndef HATKERNEL_HPP_
#define HATKERNEL_HPP_

#include "spipack/Tools/Kernels/CompactKernel.hpp"

namespace spi {
namespace Tools {

/// An implementation of the hat kernel.
/**
  The hat kernel is a compact kernel (spi::Tools::CompactKernel) such that
  \f{equation*}{
    k(\theta) = \begin{cases}
      1 & \mbox{if } \theta \in [0,1] \\
      0 & \mbox{else.}
    \end{cases}
  \f}
*/
class HatKernel : public CompactKernel {
public:
  /// Construct a hat kernel
  HatKernel();

  virtual ~HatKernel() = default;
protected:

  /// Evaluate the hat kernel function \f$k(\theta)\f$
  /**
    We have already checked that \f$\theta \leq 1\f$---this function implements the support of the kernel. Therefore, we just need to return \f$1\f$.
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateCompactKernel(double const theta) const override;

private:
};

} // namespace Tools
} // namespace spi

#endif
