#include "ComponentwiseQuadraticExternalAcceleration.hpp"

ComponentwiseQuadraticExternalAcceleration::ComponentwiseQuadraticExternalAcceleration(YAML::Node const& options) : VaryExternalAcceleration(options) {}

Eigen::VectorXd ComponentwiseQuadraticExternalAcceleration::ExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& vel, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const {
  // return |0-v| * (0-v) = -|v| * v
  return -vel.array().abs()*vel.array();
}
