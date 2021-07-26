// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "PlotQuantitiesOfInterest.h"

#include "spipack/KineticEquations/KineticModels.hpp"

using namespace spi::KineticEquations;

spiEX_Boltzmann::PlotQuantitiesOfInterest::PlotQuantitiesOfInterest(spiEX_Boltzmann::BoltzmannSolver& solver) {
  // @TODO Please insert your code here.
}

spiEX_Boltzmann::PlotQuantitiesOfInterest::~PlotQuantitiesOfInterest() {
}

void spiEX_Boltzmann::PlotQuantitiesOfInterest::startPlotting( double time) {}

void spiEX_Boltzmann::PlotQuantitiesOfInterest::finishPlotting() {}

void spiEX_Boltzmann::PlotQuantitiesOfInterest::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  // extract the mass density
  const double massDensity = std::max(1.0e-10, Q[0]);

  // extract the velocity
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> velocity;
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    velocity(i) = Q[i+1]/massDensity;
  }

  // extract the energy
  const double energy = std::max(0.0, Q[3])/massDensity;

  // extract the location
  const Eigen::Map<const Eigen::Matrix<double , MacroscaleInformation::dim, 1> > loc(&x[0], MacroscaleInformation::dim);

  // get the kinetic models
  auto kineticModels = KineticModels::ExaHyPEKineticModels();
  assert(kineticModels);

  // compute the flux matrix
  const Eigen::Matrix<double, MacroscaleInformation::dim, 4> flux = kineticModels->FluxMatrix(loc, massDensity, velocity, energy);

  outputQuantities[0] = flux(0, 1);
  outputQuantities[1] = flux(0, 2);
  outputQuantities[2] = flux(1, 2);

  const int writtenUnknowns = 6;
  for (int i=3; i<writtenUnknowns; i++){
    outputQuantities[i] = (double)i;
  }
}
