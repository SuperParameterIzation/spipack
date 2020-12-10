#include "QuadraticExternalAcceleration.hpp"

QuadraticExternalAcceleration::QuadraticExternalAcceleration(YAML::Node const& options) : VaryExternalAcceleration(options) {}

Eigen::VectorXd QuadraticExternalAcceleration::ExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& vel, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const {
  // return ||0-v|| (0-v) = -||v|| v
  return -vel.norm()*vel;
}
