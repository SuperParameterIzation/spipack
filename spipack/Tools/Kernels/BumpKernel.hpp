#ifndef BUMPKERNEL_HPP_
#define BUMPKERNEL_HPP_

#include "spipack/Tools/Kernels/CompactKernel.hpp"

namespace spi {
namespace Tools {

/// In implementation of the bump kernel
/**
  The bump kernel is a compact kernel (spi::Tools::CompactKernel) such that
  \f{equation*}{
    k(\theta) = \begin{cases}
      a_0 \exp{\left( a_1 \left( 1 - \frac{1}{1-\theta^p} \right) \right) } & \mbox{if } \theta \in [0,1) \\
      0 & \mbox{else}
    \end{cases}
  \f}
  given the parameters \f$a_0, a_1>0\f$ and the exponent \f$p>0\f$.

  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Magnitude"   | <tt>double</tt> | <tt>1.0</tt> | The magnitude of the kernel (the parameter \f$a_0\f$). |
      "Scale"   | <tt>double</tt> | <tt>1.0</tt> | The scale of the kernel (the parameter \f$a_1\f$). |
      "Exponent"   | <tt>double</tt> | <tt>1.0</tt> | The exponent of the kernel (the parameter \f$p\f$). |
*/
class BumpKernel : public CompactKernel {
public:

  /// Construct a bump kernel
  /**
  @param[in] options Options for this kernel function
  */
  BumpKernel(YAML::Node const& options);

  virtual ~BumpKernel() = default;

protected:

  /// Evaluate the hat kernel function \f$k(\theta)\f$
  /**
    @param[in] theta The value of \f$\theta\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateCompactKernelImpl(double const theta) const override;

private:

  /// The magnitude of the bump kernel (the parameter \f$a_0\f$)
  /**
  Defaults to \f$1\f$.
  */
  const double mag;

  /// The scale of the bump kernel (the parameter \f$a_1\f$)
  /**
  Defaults to \f$1\f$.
  */
  const double scale;

  /// The exponent of the bump kernel (the parameter \f$p\f$)
  /**
  Defaults to \f$1\f$.
  */
  const double expon;
};

} // namespace Tools
} // namespace spi

#endif
