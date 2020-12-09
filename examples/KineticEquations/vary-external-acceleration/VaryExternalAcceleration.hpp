#ifndef VARYEXTERNALACCELERATION_HPP_
#define VARYEXTERNALACCELERATION_HPP_

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/KineticEquations/ConditionalVelocityDistribution.hpp>

class VaryExternalAcceleration : public spi::KineticEquations::ConditionalVelocityDistribution {
public:
  /// The dimension of the state space
  inline static constexpr std::size_t dim = 2;

  /// Construct the conditional velocity distribution with a linear external forcing
  /**
  @param[in] options Setup options
  */
  inline VaryExternalAcceleration(YAML::Node const& options) : ConditionalVelocityDistribution(
    Eigen::Matrix<double, dim, 1>::Zero(dim),
    std::make_shared<muq::Modeling::Gaussian>(dim)->AsVariable(), DefaultMacroscaleInformation(),
    options)
    {}

  virtual ~VaryExternalAcceleration() = default;

  /**
  \return The macro-scale information with constat mass density, zero velocity, and zero derivatives.
  */
  inline static std::shared_ptr<MacroscaleInformation> DefaultMacroscaleInformation() {
      return std::make_shared<MacroscaleInformation>(1.0, Eigen::Matrix<double, dim, 1>::Zero(dim, 1), 0.0, Eigen::Matrix<double, dim, 1>::Zero(dim, 1));
    }
private:
};

#endif
