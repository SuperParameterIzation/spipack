#ifndef KERNEL_HPP_
#define KERNEL_HPP_

#include <yaml-cpp/yaml.h>

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

  /// Evaluate the kernel function
  /**
    @param[in] theta The (non-negative) argument to the kernel function
  */
  //virtual double operator()(double const theta) const = 0;
private:
};

} // namespace NumericalSolvers
} // namespace spi

#endif
