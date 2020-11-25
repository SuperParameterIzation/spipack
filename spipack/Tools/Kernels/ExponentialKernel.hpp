#ifndef EXPONENTIALKERNEL_HPP_
#define EXPONENTIALKERNEL_HPP_

#include "spipack/Tools/Kernels/IsotropicKernel.hpp"

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
      "Magnitude"   | <tt>double</tt> | <tt>1.0</tt> | The magnitude of the kernel (the parameter \f$a_0\f$). |
      "Scale"   | <tt>double</tt> | <tt>0.5</tt> | The parameter \f$a_1\f$. |
      "Exponent"   | <tt>double</tt> | <tt>1.0</tt> | The value of the exponent (the parameter \f$p\f$). |
*/
class ExponentialKernel : public IsotropicKernel {
public:

  /// Construct an exponential kernel
  /**
    @param[in] options Options for this kernel function
  */
  ExponentialKernel(YAML::Node const& options);

  virtual ~ExponentialKernel() = default;

  /// Evaluate the exponential kernel function \f$k(\theta)\f$
  /**
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateIsotropicKernel(double const theta) const override;

  /// Evaluate the derivative of the kernel function \f$\frac{d k}{d \theta}\f$
  /**
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$
    \return The kernel derivative \f$\frac{d k}{d \theta}\f$
  */
  virtual double IsotropicKernelDerivative(double const theta) const override;

  /// Evaluate the derivative of the kernel function \f$\frac{d k}{d \theta}\f$
  /**
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$
    \return The kernel derivative \f$\frac{d k}{d \theta}\f$
  */
  virtual double IsotropicKernelSecondDerivative(double const theta) const override;

  /// Get the kernel's magnitude parameter \f$a_0\f$
  /**
  \return The kernel's magnitude parameter \f$a_0\f$
  */
  double Magnitude() const;

  /// Get the kernel's scale parameter \f$a_1\f$
  /**
  \return The kernel's scale parameter \f$a_1\f$
  */
  double Scale() const;

  /// Get the kernel's exponent \f$p\f$
  /**
  \return The kernel's exponent \f$p\f$
  */
  double Exponent() const;

  /// Compute the integral of the kernel
  /**
  Compute the integral
  \f{equation*}{
  m_n = \int_{\mathbb{R}^d} x_1^n k(\|\boldsymbol{x}\|^2) \, d\boldsymbol{x}
  \f}
  @param[in] dim The dimension of the problem \f$d\f$
  @param[in] n The parameter \f$n\f$ that defines \f$m_n\f$
  \return The numerically computed integral
  */
  virtual double Integrate(std::size_t const dim, std::size_t const n) const override;

protected:

  /// Private default constructor for the serialization
  ExponentialKernel() = default;

private:

  /// Make cereal::access a friend for serialization
  friend class cereal::access;

  /// Serialize so we can save to archive/buffers
  /**
  @param[in] ar The archive/buffer
  */
  template<class Archive>
  inline void serialize(Archive& ar) { ar(cereal::base_class<IsotropicKernel>(this), mag, scale, expon); }

  /// The magnitude parameter \f$a_0\f$
  double mag = 1.0;

  /// The scale parameter \f$a_1\f$
  double scale = 0.5;

  /// The exponent \f$p\f$
  double expon = 1.0;
};

} // namespace Tools
} // namespace spi

#endif
