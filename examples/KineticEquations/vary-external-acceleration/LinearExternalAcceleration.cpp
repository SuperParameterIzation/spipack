#include "LinearExternalAcceleration.hpp"

LinearExternalAcceleration::LinearExternalAcceleration(YAML::Node const& options) : VaryExternalAcceleration(options) {}

Eigen::VectorXd LinearExternalAcceleration::ExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& vel, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const {
  // return 0-v = -v
  return -vel;
}
