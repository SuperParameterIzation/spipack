#include "ComponentwiseQuadraticExternalAcceleration.hpp"

using namespace spi::KineticEquations;

ComponentwiseQuadraticExternalAcceleration::ComponentwiseQuadraticExternalAcceleration(YAML::Node const& options) : VaryExternalAcceleration(options) {}

Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> ComponentwiseQuadraticExternalAcceleration::ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& vel, double const time) const {
  // return |0-v| * (0-v) = -|v| * v
  return -vel.array().abs()*vel.array();
}
