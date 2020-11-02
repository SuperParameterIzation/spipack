#ifndef ISOTROPICKERNEL_HPP_
#define ISOTROPICKERNEL_HPP_

#include "spipack/Tools/Kernels/Kernel.hpp"

namespace spi {
namespace Tools {

/// An implementation of an isotropic kernel
/**
  An isotropic kernel takes the form of a decreasing function \f$k:\mathbb{R}^{+} \mapsto \mathbb{R}^{+}\f$ such that \f$k(\theta) = k(\boldsymbol{x}_1, \boldsymbol{x}_2) = k(\|\boldsymbol{x}_1 - \boldsymbol{x}_2\|^2)\f$.

  Note: we define the kernel in terms \f$\|\boldsymbol{x}_1 - \boldsymbol{x}_2\|^2\f$ (rather than \f$\|\boldsymbol{x}_1 - \boldsymbol{x}_2\|\f$) for numerical convenience.
*/
class IsotropicKernel : public Kernel {
public:

  /// Construct an isotropic kernel
  IsotropicKernel();

  virtual ~IsotropicKernel() = default;

  /// Evaluate the kernel function \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  /**
    @param[in] x1 The first argument to the kernel function
    @param[in] x2 The second argument to the kernel function
    \return The kernel evaluation \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  */
  virtual double Evaluate(Eigen::Ref<Eigen::VectorXd> const& x1, Eigen::Ref<Eigen::VectorXd> const& x2) const override;

  /// Evaluate the kernel function \f$k(\theta)\f$
  /**
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateIsotropicKernel(double const theta) const = 0;

  /// Evaluate the kernel function \f$k(\theta)\f$
  /**
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$
    \return The kernel evaluation \f$k(\theta)\f$
  */
  double operator()(double const theta) const;

private:
};

} // namespace Tools
} // namespace spi

#endif
