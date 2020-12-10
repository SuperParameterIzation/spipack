#ifndef HATKERNEL_HPP_
#define HATKERNEL_HPP_

#include <cereal/archives/binary.hpp>

#include "spipack/Tools/Kernels/CompactKernel.hpp"

namespace spi {
namespace Tools {

/// An implementation of the hat kernel.
/**
  The hat kernel is a compact kernel (spi::Tools::CompactKernel) such that
  \f{equation*}{
    k(\theta) = \begin{cases}
      a & \mbox{if } \theta \in [0,1] \\
      0 & \mbox{else}
    \end{cases}
  \f}
  given the parameter \f$a>0\f$.

  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Magnitude"   | <tt>double</tt> | <tt>1.0</tt> | The magnitude of the kernel (the parameter \f$a\f$). |
*/
class HatKernel : public CompactKernel {
public:

  /// Construct a hat kernel
  /**
    @param[in] options Options for this kernel function
  */
  HatKernel(YAML::Node const& options = YAML::Node());

  virtual ~HatKernel() = default;

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

  /// Get the kernel's magnitude parameter \f$a\f$
  /**
  \return The kernel's magnitude parameter \f$a\f$
  */
  double Magnitude() const;

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

  /// Evaluate the hat kernel function \f$k(\theta)\f$
  /**
    We have already checked that \f$\theta \leq 1\f$---this function implements the support of the kernel. Therefore, we just need to return the magnitude \f$a\f$.
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateCompactKernelImpl(double const theta) const override;

private:

  /// Make cereal::access a friend for serialization
  friend class cereal::access;

  /// Serialize so we can save to archive/buffers
  /**
  @param[in] ar The archive/buffer
  */
  template<class Archive>
  inline void serialize(Archive& ar) { ar(cereal::base_class<CompactKernel>(this), mag); }

  /// The magnitude of the kernel (the parameter \f$a\f$).
  double mag = 1.0;
};

} // namespace Tools
} // namespace spi

#endif
