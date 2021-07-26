// This file was generated by the ExaHyPE toolkit.
// It will NOT be regenerated or overwritten.
// Please adapt it to your own needs.
//
// ========================
//   www.exahype.eu
// ========================

#include "BoltzmannSolver_ADERDG.h"

#include "BoltzmannSolver_ADERDG_Variables.h"

#include "kernels/KernelUtils.h"
#include "peano/utils/Loop.h"

#include "spipack/KineticEquations/KineticModels.hpp"

using namespace spi::KineticEquations;

tarch::logging::Log spiEX_Boltzmann::BoltzmannSolver_ADERDG::_log("spiEX_Boltzmann::BoltzmannSolver_ADERDG" );

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::init(const std::vector<std::string>& cmdlineargs, const exahype::parser::ParserView& constants) {}

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::adjustPointSolution(const double* const x, const double t, const double dt, double* const Q) {
  //assert(false);
  if( tarch::la::equals(t, 0.0) ) {
    // get the kinetic models
    auto kineticModels = KineticModels::ExaHyPEKineticModels();
    assert(kineticModels);

    // the spatial location
    Eigen::Map<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > loc(x, MacroscaleInformation::dim);

    // initial mass density
    Q[0] = kineticModels->InitialMassDensity(loc);

    // initial velocity
    Eigen::Map<Eigen::Matrix<double, MacroscaleInformation::dim, 1> > vel(&Q[1], MacroscaleInformation::dim); // inital mass density gradient---only required to pass to this function
    vel = kineticModels->InitialExpectedVelocity(loc);

    // we actually need the initial momentum
    for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { vel(i) *= Q[0]; }

    // the initial expected energy
    Q[3] = kineticModels->InitialExpectedEnergy(loc)*Q[0];

    // the coordinates
    Q[4] = x[0]; Q[5] = x[1];
  }

  // check the mass and energy (they must be positive)
  Q[0] = std::max(0.0, Q[0]);
  Q[3] = std::max(0.0, Q[3]);

  // check the coordinates
  assert(Q[4]>-1.0e-10); assert(Q[4]<1.0+1.0e-10);
  assert(Q[5]>-1.0e-10); assert(Q[5]<1.0+1.0e-10);
}

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::boundaryValues(
  const double* const x,
  const double t,
  const double dt,
  const int faceIndex,
  const int direction,
  const double* const fluxIn,
  const double* const stateIn,
  const double* const gradStateIn,
  double* const fluxOut,
  double* const stateOut)
{
  assert(false);
  assert(stateIn[4]>-1.0e-3); assert(stateIn[4]<1.0+1.0e-3);
  assert(stateIn[5]>-1.0e-3); assert(stateIn[5]<1.0+1.0e-3);

  // copy the state
  std::copy_n(stateIn, NumberOfVariables+NumberOfParameters, stateOut);

  // for reflecting boundaries
  stateOut[1+direction] = -stateOut[1+direction];

  // no-flow boundary conditions
  //stateOut[1] = 0.0;
  //stateOut[2] = 0.0;

  // compute the flux using stateOut
  double _F[DIMENSIONS][NumberOfVariables]={0.0};
  #if DIMENSIONS==2
  double* F[2] = {_F[0], _F[1]};
  #elif DIMENSIONS==3
  double* F[3] = {_F[0], _F[1], _F[2]};
  #endif
  flux(stateOut, F);
  std::copy_n(F[direction], NumberOfVariables, fluxOut);
}

exahype::solvers::Solver::RefinementControl spiEX_Boltzmann::BoltzmannSolver_ADERDG::refinementCriterion(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& cellCentre, const tarch::la::Vector<DIMENSIONS,double>& cellSize, double t, const int level) { return exahype::solvers::Solver::RefinementControl::Keep; }

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::eigenvalues(const double* const Q, const int direction, double* const lambda) {
  // get the kinetic models
  auto kineticModels = KineticModels::ExaHyPEKineticModels();
  assert(kineticModels);

  // extract the location
  const Eigen::Map<const Eigen::Matrix<double , MacroscaleInformation::dim, 1> > loc(&Q[4], MacroscaleInformation::dim);

  assert(Q[4]>-1.0e-3); assert(Q[4]<1.0+1.0e-3);
  assert(loc(0)>-1.0e-3); assert(loc(0)<1.0+1.0e-3);
  assert(Q[5]>-1.0e-3); assert(Q[5]<1.0+1.0e-3);
  assert(loc(1)>-1.0e-3); assert(loc(1)<1.0+1.0e-3);

  // extract the mass density
  const KineticModels::ADScalar massDensity(2+MacroscaleInformation::dim, 0, std::max(1.0e-10, Q[0]));

  // extract the velocity
  Eigen::Matrix<KineticModels::ADScalar, MacroscaleInformation::dim, 1> velocity;
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { velocity(i) = KineticModels::ADScalar(2+MacroscaleInformation::dim, i+1, Q[i+1])/massDensity; }

  // extract the energy
  const KineticModels::ADScalar energy = KineticModels::ADScalar(2+MacroscaleInformation::dim, 1+MacroscaleInformation::dim, std::max(0.0, Q[3]))/massDensity;

  // compute the flux matrix
  const Eigen::Matrix<KineticModels::ADScalar, MacroscaleInformation::dim, 4> flux = kineticModels->FluxMatrix(loc, massDensity, velocity, energy);
  assert(flux.cols()==NumberOfVariables);

  Eigen::Matrix<double, NumberOfVariables, NumberOfVariables> dFdU = Eigen::Matrix<double, NumberOfVariables, NumberOfVariables>::Zero(NumberOfVariables, NumberOfVariables);

  for( std::size_t i=0; i<NumberOfVariables; ++i ) {
    for( std::size_t j=0; j<NumberOfVariables; ++j ) { dFdU(i, j) = flux(direction, j).dx(i); }
  }

  Eigen::Map<Eigen::Matrix<double, NumberOfVariables, 1> > eigs(lambda, NumberOfVariables);
  eigs = dFdU.eigenvalues().real();
  if( eigs.array().abs().maxCoeff()<1.0e-10 ) { eigs = dFdU.transpose().eigenvalues().real(); }
}

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::flux(const double* const Q, double** const F) {
  assert(false);
  // get the kinetic models
  auto kineticModels = KineticModels::ExaHyPEKineticModels();
  assert(kineticModels);

  // extract the location
  const Eigen::Map<const Eigen::Matrix<double , MacroscaleInformation::dim, 1> > loc(&Q[4], MacroscaleInformation::dim);

  assert(Q[4]>-1.0e-3); assert(Q[4]<1.0+1.0e-3);
  assert(loc(0)>-1.0e-3); assert(loc(0)<1.0+1.0e-3);
  assert(Q[5]>-1.0e-3); assert(Q[5]<1.0+1.0e-3);
  assert(loc(1)>-1.0e-3); assert(loc(1)<1.0+1.0e-3);

  // extract the mass density
  const double massDensity = std::max(1.0e-10, Q[0]);

  // extract the velocity
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> velocity;
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { velocity(i) = Q[i+1]/massDensity; }

  // extract the energy
  const double energy = std::max(0.0, Q[3])/massDensity;

  // compute the flux matrix
  const Eigen::Matrix<double, MacroscaleInformation::dim, 4> flux = kineticModels->FluxMatrix(loc, massDensity, velocity, energy);
  assert(flux.cols()==NumberOfVariables);

  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    for( std::size_t j=0; j<NumberOfVariables; ++j ) { F[i][j] = flux(i, j); }
  }
}

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  assert(false);
  //for( std::size_t i=0; i<NumberOfVariables; ++i ) { assert(!std::isnan(Q[i])); }
  assert(std::abs(x[0]-Q[4])<1.0e-10);
  assert(std::abs(x[1]-Q[5])<1.0e-10);

  const double mu = std::max(1.0e-10, Q[0]);

  const double f1 = std::abs(1.0-Q[1]/mu)*(1.0-Q[1]/mu);
  const double f2 = std::abs(0.0-Q[2]/mu)*(0.0-Q[2]/mu);

  const double p = Q[3]-0.5*(Q[1]*Q[1] + Q[2]*Q[2])/(mu*mu);
  //const double p = Q[3];
  const double k = 25.0;

  //std::cout << "p: " << p << " " << 1.0/(1.0+std::exp(-2.0*k*(p-1.25))) << std::endl;

  S[0] = 0.0;
  S[1] = Q[0]*f1;
  S[2] = Q[0]*f2;
  S[3] = (f1*Q[1] + f2*Q[2])/mu;
  //S[3] += 1.0/(1.0+std::exp(-2.0*k*(p-0.5)));
  //S[3] += p*std::exp(-0.5*((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5)));
}

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::nonConservativeProduct(const double* const Q, const double* const gradQ, double* const BgradQ) {
  assert(false);
  /*for( std::size_t i=0; i<10; ++i ) {
    std::cout << "(ADERDG) i: " << i << " gradQ[i]: " << gradQ[i] << std::endl;
  }
  std::cout << std::endl;*/

  kernels::idx2 idx_gradQ(DIMENSIONS, NumberOfVariables+NumberOfParameters);

  const double mu = std::max(1.0e-10, Q[0]);
  const double e = std::max(0.0, Q[3])/mu;
  const double mugradnorm = std::sqrt(gradQ[idx_gradQ(0, 0)]*gradQ[idx_gradQ(0, 0)]+gradQ[idx_gradQ(1, 0)]*gradQ[idx_gradQ(1, 0)]);
  const double p = e-0.5*(Q[1]*Q[1] + Q[2]*Q[2])/(mu*mu);

  const double divMomentum = gradQ[idx_gradQ(0, 1)] + gradQ[idx_gradQ(1, 2)];

  assert(false);

  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  //BgradQ[3] = e*mugradnorm*(1.0 + (gradQ[idx_gradQ(0, 1)] + gradQ[idx_gradQ(1, 2)] - Q[1]*gradQ[idx_gradQ(0, 0)]/Q[0]  -Q[2]*gradQ[idx_gradQ(1, 0)]/Q[0]));
  //BgradQ[3] = 3.0*(gradQ[idx_gradQ(0, 0)]+gradQ[idx_gradQ(1, 0)]);
  //BgradQ[3] = mugradnorm*e;

  //BgradQ[3] = e*mu*divMomentum;

}

void spiEX_Boltzmann::BoltzmannSolver_ADERDG::mapDiscreteMaximumPrincipleObservables(double* const observables, const double* const Q) const {
  for (int i=0; i<NumberOfVariables; ++i) {
    observables[i] = Q[i];
  }
}

bool spiEX_Boltzmann::BoltzmannSolver_ADERDG::isPhysicallyAdmissible(
        const double* const Q,
        const double* const observablesMin,const double* const observablesMax,
        const bool wasTroubledInPreviousTimeStep,
        const tarch::la::Vector<DIMENSIONS,double>& center,
        const tarch::la::Vector<DIMENSIONS,double>& dx,
        const double t) const {
          return false;
        }
