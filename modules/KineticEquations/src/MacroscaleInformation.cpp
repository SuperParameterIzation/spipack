#include "spipack/KineticEquations/MacroscaleInformation.hpp"

using namespace spi::KineticEquations;

MacroscaleInformation::MacroscaleInformation() {}

MacroscaleInformation::MacroscaleInformation(double const massDensity, Eigen::Matrix<double, dim, 1> const& logMassDensityGrad, Eigen::Matrix<double, dim, 1> const& velocity, double const velocityDiv, Coordinates const& coordinates) :
massDensity(massDensity),
velocity(velocity),
logMassDensityGrad(logMassDensityGrad),
velocityDiv(velocityDiv),
coordinates(coordinates)
{}
