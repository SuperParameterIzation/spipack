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
      a & \mbox{if } \theta \in [0,1] \\
      0 & \mbox{else.}
    \end{cases}
  \f}
  given the parameter \f$a>0\f$.

  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Magnitude"   | double | <tt>1.0</tt> | The magnitude of the kernel (the parameter \f$a\f$). |
*/
class HatKernel : public CompactKernel {
public:

  /// Construct a hat kernel
  /**
    @param[in] options Options for this kernel function
  */
  HatKernel(YAML::Node const& options);

  virtual ~HatKernel() = default;

  /// Evaluate the hat kernel function \f$k(\theta)\f$
  /**
    We have already checked that \f$\theta \leq 1\f$---this function implements the support of the kernel. Therefore, we just need to return the magnitude \f$a\f$.
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateCompactKernel(double const theta) const override;

private:

  /// The magnitude of the kernel (the parameter \f$a\f$).
  const double mag;
};

} // namespace Tools
} // namespace spi

#endif
