#include "spipack/Tools/Kernels/Kernel.hpp"

using namespace spi::Tools;

Kernel::Kernel() {}

double Kernel::operator()(Eigen::Ref<Eigen::VectorXd> const& x1, Eigen::Ref<Eigen::VectorXd> const& x2) const { return Evaluate(x1, x2); }
