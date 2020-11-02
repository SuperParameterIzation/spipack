#ifndef COMPACTKERNEL_HPP_
#define COMPACTKERNEL_HPP_

#include "spipack/Tools/Kernels/IsotropicKernel.hpp"

namespace spi {
namespace Tools {

/// An implementation of a compact kernel
/**
  A compact kernel takes the form of a decreasing function \f$k:\mathbb{R}^{+} \mapsto \mathbb{R}^{+}\f$ such that \f$k(\theta) = 0\f$ if \f$\theta \notin [0,1]\f$ (i.e., the support of the kernel is \f$[0,1]\f$.)
*/
class CompactKernel : public IsotropicKernel {
public:

  /// Construct a compact kernel
  CompactKernel();

  virtual ~CompactKernel() = default;

  /// Evaluate the kernel function \f$k(\theta)\f$
  /**
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateIsotropicKernel(double const theta) const override;

protected:

  /// Evaluate the kernel function \f$k(\theta)\f$
  /**
    We have already checked that \f$\theta \leq 1\f$---this function implements the support of the kernel.
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateCompactKernel(double const theta) const = 0;

private:
};


} // namespace Tools
} // namespace spi

#endif
