#include "spipack/Tools/Kernels/HatKernel.hpp"

using namespace spi::Tools;

HatKernel::HatKernel() : CompactKernel() {}

double HatKernel::EvaluateCompactKernel(double const theta) const { return 1.0; }
