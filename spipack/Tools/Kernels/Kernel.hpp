#ifndef KERNEL_HPP_
#define KERNEL_HPP_

#include <Eigen/Core>

namespace spi {
namespace Tools {

/// The kernel function implementation
/**
  This class implements a kernel function \f$k: \mathbb{R}^{d} \times \mathbb{R}^{d} \mapsto \mathbb{R}^{+}\f$.
*/
class Kernel {
public:

  /// Construct the kernel function
  Kernel();

  virtual ~Kernel() = default;

  /// Evaluate the kernel function \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  /**
    @param[in] x1 The first argument to the kernel function
    @param[in] x2 The second argument to the kernel function
    \return The kernel evaluation \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  */
  virtual double Evaluate(Eigen::Ref<Eigen::VectorXd> const& x1, Eigen::Ref<Eigen::VectorXd> const& x2) const = 0;

  /// Evaluate the kernel function  \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  /**
    @param[in] x1 The first argument to the kernel function
    @param[in] x2 The second argument to the kernel function
    \return The kernel evaluation \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  */
  double operator()(Eigen::Ref<Eigen::VectorXd> const& x1, Eigen::Ref<Eigen::VectorXd> const& x2) const;

private:
};

} // namespace Tools
} // namespace spi

#endif
