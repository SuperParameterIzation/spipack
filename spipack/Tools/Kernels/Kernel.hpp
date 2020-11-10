#ifndef KERNEL_HPP_
#define KERNEL_HPP_

#include <iostream>
#include <map>

#include <boost/function.hpp>

#include <yaml-cpp/yaml.h>

#include <Eigen/Core>

#include <MUQ/Utilities/RegisterClassName.h>

namespace spi {
namespace Tools {

/// The kernel function implementation
/**
  This class implements a kernel function \f$k: \mathbb{R}^{d} \times \mathbb{R}^{d} \mapsto \mathbb{R}^{+}\f$.
*/
class Kernel {
public:

  /// Construct the kernel function
  /**
  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "Kernel"   | <tt>std::string</tt> | - | The name of the kernel function   |
  */
  Kernel(YAML::Node const& options);

  virtual ~Kernel() = default;

  /// A static constructor to create a kernel function
  /**
    @param[in] options Options for this constructor
  */
  static std::shared_ptr<Kernel> Construct(YAML::Node const& options);

  /// The constructor type to create a kernel function
  typedef boost::function<std::shared_ptr<Kernel>(YAML::Node)> KernelConstructor;

  /// A map from the kernel name to its corresponding constructor
  typedef std::map<std::string, KernelConstructor> ConstructKernelMap;

  /// Get the map from kernel name to its constructor
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
  virtual double Evaluate(Eigen::Ref<const Eigen::VectorXd> const& x1, Eigen::Ref<const Eigen::VectorXd> const& x2) const = 0;

  /// Evaluate the kernel function  \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  /**
    @param[in] x1 The first argument to the kernel function
    @param[in] x2 The second argument to the kernel function
    \return The kernel evaluation \f$k(\boldsymbol{x}_1, \boldsymbol{x}_2)\f$
  */
  double operator()(Eigen::Ref<const Eigen::VectorXd> const& x1, Eigen::Ref<const Eigen::VectorXd> const& x2) const;

private:
};

} // namespace Tools
} // namespace spi

#define SPIPACK_REGISTER_KERNEL(NAME) static auto reg ##NAME		\
  = spi::Tools::Kernel::KernelMap()->insert(std::make_pair(#NAME, muq::Utilities::shared_factory<NAME>()));

#endif
