#include "LinearExternalAcceleration.hpp"

using namespace spi::KineticEquations;

LinearExternalAcceleration::LinearExternalAcceleration(YAML::Node const& options) : VaryExternalAcceleration(options) {}

Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> LinearExternalAcceleration::ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& vel, double const time) const {
  // return 0-v = -v
  return -vel;
}
