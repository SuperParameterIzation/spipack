#ifndef LINEAREXTERNALACCELERATION_HPP_
#define LINEAREXTERNALACCELERATION_HPP_

#include "VaryExternalAcceleration.hpp"

class LinearExternalAcceleration : public VaryExternalAcceleration {
public:
  /// Construct the conditional velocity distribution with a linear external forcing
  /**
  @param[in] options Setup options
  */
  LinearExternalAcceleration(YAML::Node const& options);

  virtual ~LinearExternalAcceleration() = default;

  /**
  @param[in] vel The particle velocity \$\boldsymbol{v}\f$
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  virtual Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& vel, double const time) const override;

private:
};

#endif
