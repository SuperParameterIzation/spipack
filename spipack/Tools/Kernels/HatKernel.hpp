#ifndef HATKERNEL_HPP_
#define HATKERNEL_HPP_

#include "spipack/Tools/Kernels/Kernel.hpp"

namespace spi {
namespace NumericalSolvers {

/// An implementation of the hat kernel for the graph Laplacian
class HatKernel {
  virtual ~HatKernel() = default;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
