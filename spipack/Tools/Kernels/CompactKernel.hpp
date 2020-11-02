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
  /**
    @param[in] options Options for this kernel function
  */
  CompactKernel(YAML::Node const& options);

  virtual ~CompactKernel() = default;

  // A static constructor to create a compact kernel function
  /**
    @param[in] options Options for this constructor
  */
  static std::shared_ptr<CompactKernel> Construct(YAML::Node const& options);

  /// The constructor type to create a compact kernel function
  typedef boost::function<std::shared_ptr<CompactKernel>(YAML::Node)> KernelConstructor;

  /// A map from the compact kernel name to its corresponding constructor
  typedef std::map<std::string, KernelConstructor> ConstructKernelMap;

  /// Get the map from isotropic kernel name to its constructor
  /**
    \return The map from kernel name to its constructor
  */
  static std::shared_ptr<ConstructKernelMap> KernelMap();

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

#define SPIPACK_REGISTER_COMPACT_KERNEL(NAME) static auto regCompact ##NAME		\
  = spi::Tools::CompactKernel::KernelMap()->insert(std::make_pair(#NAME, muq::Utilities::shared_factory<NAME>()));

} // namespace Tools
} // namespace spi

#endif
