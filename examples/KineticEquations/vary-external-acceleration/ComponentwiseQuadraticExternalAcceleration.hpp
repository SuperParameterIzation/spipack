#ifndef COMPONENTWISEQUADRATICEXTERNALACCELERATION_HPP_
#define COMPONENTWISEQUADRATICEXTERNALACCELERATION_HPP_

#include "VaryExternalAcceleration.hpp"

class ComponentwiseQuadraticExternalAcceleration : public VaryExternalAcceleration {
public:
  /// Construct the conditional velocity distribution with a componentwise quadratic external forcing
  /**
  @param[in] options Setup options
  */
  ComponentwiseQuadraticExternalAcceleration(YAML::Node const& options);

  virtual ~ComponentwiseQuadraticExternalAcceleration() = default;

  /**
  @param[in] vel The particle velocity \$\boldsymbol{v}\f$
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  virtual Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& vel, double const time) const override;

private:
};

#endif
