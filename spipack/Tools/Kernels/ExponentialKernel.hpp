#ifndef EXPONENTIALKERNEL_HPP_
#define EXPONENTIALKERNEL_HPP_

#include "spipack/Tools/Kernels/CompactKernel.hpp"

namespace spi {
namespace Tools {

/// An implementation of the exponential kernel.
/**
  The hat kernel is an isotropic kernel (spi::Tools::IsotropicKernel) such that
  \f{equation*}{
    k(\theta) = a_0 \exp{(-a_1 \vert \theta \vert^{2p})}
  \f}
  given parameters \f$a_0, a_1>0\f$ and exponent \f$p>0\f$.

  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Magnitude"   | double | <tt>1.0</tt> | The magnitude of the kernel (the parameter \f$a_0\f$). |
      "Scale"   | double | <tt>1.0</tt> | The parameter \f$a_1\f$. |
      "Exponent"   | double | <tt>1.0</tt> | The value of the exponent (the parameter \f$p\f$). |
*/
class ExponentialKernel : public IsotropicKernel {
public:

  /// Construct an exponential kernel
  /**
    @param[in] options Options for this kernel function
  */
  ExponentialKernel(YAML::Node const& options);

  virtual ~ExponentialKernel() = default;

private:

  /// Evaluate the hat kernel function \f$k(\theta)\f$
  /**
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateIsotropicKernel(double const theta) const override;

  /// The magnitude parameter \f$a_0\f$
  const double mag;

  /// The scale parameter \f$a_1\f$
  const double scale;

  /// The exponent \f$p\f$
  const double expon;
};

} // namespace Tools
} // namespace spi

#endif
