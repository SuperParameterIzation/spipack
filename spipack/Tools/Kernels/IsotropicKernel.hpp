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
  /**
    @param[in] options Options for this kernel function
  */
  IsotropicKernel(YAML::Node const& options);

  virtual ~IsotropicKernel() = default;

  /// A static constructor to create an isotropic kernel function
  /**
    @param[in] options Options for this constructor
  */
  static std::shared_ptr<IsotropicKernel> Construct(YAML::Node const& options);

  /// The constructor type to create an isotropic kernel function
  typedef boost::function<std::shared_ptr<IsotropicKernel>(YAML::Node)> KernelConstructor;

  /// A map from the isotropic kernel name to its corresponding constructor
  typedef std::map<std::string, KernelConstructor> ConstructKernelMap;

  /// Get the map from isotropic kernel name to its constructor
  /**
    \return The map from kernel name to its constructor
  */
  static std::shared_ptr<ConstructKernelMap> KernelMap();

  /// Evaluate the kernel function \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  /**
    @param[in] x1 The first argument to the kernel function
    @param[in] x2 The second argument to the kernel function
    \return The kernel evaluation \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  */
  virtual double Evaluate(Eigen::Ref<const Eigen::VectorXd> const& x1, Eigen::Ref<const Eigen::VectorXd> const& x2) const override;

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

#define SPIPACK_REGISTER_ISOTROPIC_KERNEL(NAME) static auto regIsotropic ##NAME		\
  = spi::Tools::IsotropicKernel::KernelMap()->insert(std::make_pair(#NAME, muq::Utilities::shared_factory<NAME>()));

#endif
