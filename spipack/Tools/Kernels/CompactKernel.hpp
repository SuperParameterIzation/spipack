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

  /// Evaluate the kernel function \f$k(\theta)\f$
  /**
    We have already checked that \f$\theta \leq 1\f$---this function implements the support of the kernel.
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  double EvaluateCompactKernel(double const theta) const;

  /// Numerically compute the integral of the kernel
  /**
  Compute the integral
  \f{equation*}{
  m_n = \int_{\mathbb{R}^d} x_1^n k(\|\boldsymbol{x}\|^2) \, d\boldsymbol{x}
  \f}
  This default implementation does this numerically using Gauss-Legendre quadrature and a full tensor product.
  @param[in] dim The dimension of the problem \f$d\f$
  @param[in] n The parameter \f$n\f$ that defines \f$m_n\f$
  \return The numerically computed integral
  */
  virtual double NumericallyIntegrate(std::size_t const dim, std::size_t const n) const override;

protected:

  /// Evaluate the kernel function \f$k(\theta)\f$
  /**
    We have already checked that \f$\theta \leq 1\f$---this function implements the support of the kernel.
    @param[in] theta The value of \f$\theta = \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2\f$ (note that \f$0 \leq \theta \leq 1\f$)
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double EvaluateCompactKernelImpl(double const theta) const = 0;

  /// Private default constructor for the serialization
  CompactKernel() = default;

private:

  /// Make cereal::access a friend for serialization
  friend class cereal::access;

  /// Serialize so we can save to archive/buffers
  /**
  @param[in] ar The archive/buffer
  */
  template<class Archive>
  inline void serialize(Archive& ar) { ar(cereal::base_class<IsotropicKernel>(this)); }
};

#define SPIPACK_REGISTER_COMPACT_KERNEL(NAME) static auto regCompact ##NAME		\
  = spi::Tools::CompactKernel::KernelMap()->insert(std::make_pair(#NAME, muq::Utilities::shared_factory<NAME>()));

} // namespace Tools
} // namespace spi

#endif
