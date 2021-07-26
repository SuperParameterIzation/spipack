#include "QuadraticExternalAcceleration.hpp"

using namespace spi::KineticEquations;

QuadraticExternalAcceleration::QuadraticExternalAcceleration(YAML::Node const& options) : VaryExternalAcceleration(options) {}

Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> QuadraticExternalAcceleration::ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& vel, double const time) const {
  // return ||0-v|| (0-v) = -||v|| v
  return -vel.norm()*vel;
}
