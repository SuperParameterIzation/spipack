#include "BoltzmannSolver_FV.h"

#include "BoltzmannSolver_FV_Variables.h"

#include "kernels/KernelUtils.h"

#include "spipack/KineticEquations/KineticModels.hpp"

using namespace spi::KineticEquations;

tarch::logging::Log spiEX_Boltzmann::BoltzmannSolver_FV::_log("spiEX_Boltzmann::BoltzmannSolver_FV" );

void spiEX_Boltzmann::BoltzmannSolver_FV::init(const std::vector<std::string>& cmdlineargs, const exahype::parser::ParserView& constants) {}

void spiEX_Boltzmann::BoltzmannSolver_FV::adjustSolution(const double* const x, const double t,const double dt, double* const Q) {
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

void spiEX_Boltzmann::BoltzmannSolver_FV::eigenvalues(const double* const Q, const int direction, double* const lambda) {
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
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    velocity(i) = KineticModels::ADScalar(2+MacroscaleInformation::dim, i+1, Q[i+1])/massDensity;
  }

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

void spiEX_Boltzmann::BoltzmannSolver_FV::boundaryValues(
  const double* const x,
  const double t,
  const double dt,
  const int faceIndex,
  const int direction,
  const double* const stateIn,
  double* const stateOut)
{
  assert(stateIn[4]>-1.0e-3); assert(stateIn[4]<1.0+1.0e-3);
  assert(stateIn[5]>-1.0e-3); assert(stateIn[5]<1.0+1.0e-3);

  // copy the state
  std::copy_n(stateIn, NumberOfVariables+NumberOfParameters, stateOut);

  // for reflecting boundaries
  stateOut[1+direction] = -stateOut[1+direction];

  // no-flow boundary conditions
  //stateOut[1] = 0.0;
  //stateOut[2] = 0.0;
}

void spiEX_Boltzmann::BoltzmannSolver_FV::flux(const double* const Q, double** const F) {
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

void spiEX_Boltzmann::BoltzmannSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // get the kinetic models
  auto kineticModels = KineticModels::ExaHyPEKineticModels();
  assert(kineticModels);

  // extract the location
  const Eigen::Map<const Eigen::Matrix<double , MacroscaleInformation::dim, 1> > loc(&Q[4], MacroscaleInformation::dim);

  assert(Q[4]>-1.0e-3); assert(Q[4]<1.0+1.0e-3);
  assert(loc(0)>-1.0e-3); assert(loc(0)<1.0+1.0e-3);
  assert(Q[5]>-1.0e-3); assert(Q[5]<1.0+1.0e-3);
  assert(loc(1)>-1.0e-3); assert(loc(1)<1.0+1.0e-3);

  //for( std::size_t i=0; i<NumberOfVariables; ++i ) { assert(!std::isnan(Q[i])); }
  //assert(std::abs(x[0]-Q[4])<1.0e-10);
  //assert(std::abs(x[1]-Q[5])<1.0e-10);

  // extract the mass density
  const double mu = std::max(1.0e-10, Q[0]);

  // extract the velocity
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> velocity;
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { velocity(i) = Q[i+1]/mu; }

  Eigen::Matrix<double, DIMENSIONS, 1> acceleration;
  double squaredDistanceToModel;
  std::tie(acceleration, squaredDistanceToModel) = kineticModels->ExpectedExternalAcceleration(loc, velocity, t);

  //Eigen::Matrix<double, DIMENSIONS, 1> accel = Eigen::Matrix<double, DIMENSIONS, 1>::Zero(DIMENSIONS, 1);
  //accel(0) = 1.0;
  //accel(0) = (loc(0)<0.5? -1.0 : 1.0);
  /*accel(0) = std::abs(accel(0)-velocity(0))*(accel(0)-velocity(0));
  accel(1) = std::abs(accel(1)-velocity(1))*(accel(1)-velocity(1));

  const double scale = kineticModels->maxModelDistanceSquared/squaredDistanceToModel;
  acceleration = scale*accel + (1.0-scale)*acceleration;*/

  //acceleration = accel;

  //Eigen::Matrix<double, DIMENSIONS, 1> externalVel;
  //externalVel(0) = 1.0; externalVel(1) = 0.0;
  //externalVel(0) = std::sin(2.0*M_PI*loc(0))*std::cos(2.0*M_PI*loc(1));
  //externalVel(1) = std::cos(2.0*M_PI*loc(0))*std::sin(2.0*M_PI*loc(1));

  //const double f1 = std::abs(externalVel(0)-Q[1]/mu)*(externalVel(0)-Q[1]/mu);
  //const double f2 = std::abs(externalVel(1)-Q[2]/mu)*(externalVel(1)-Q[2]/mu);
  //const double f1 = 1.0, f2 = 0.0;

  //const double p = Q[3]/mu-0.5*(Q[1]*Q[1] + Q[2]*Q[2])/(mu*mu);
  //const double p = Q[3];
  //const double k = 25.0;
  //const double nrgsource = -std::exp(-0.5*((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5))/10.0);

  S[0] = 0.0;
  S[1] = mu*acceleration(0);
  S[2] = mu*acceleration(1);
  S[3] = mu*acceleration.dot(velocity);
  //S[3] += 1.0/(1.0+std::exp(-2.0*k*(p-0.5)));
  //S[3] += -p;
  //S[3] += (std::sin(2.0*M_PI*t)-1.0)*p/mu*std::exp(-0.5*((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5))/1.0e-2);
}

void  spiEX_Boltzmann::BoltzmannSolver_FV::nonConservativeProduct(const double* const Q, const double* const gradQ, double* const BgradQ) {
  /*for( std::size_t i=0; i<10; ++i ) {
    std::cout << "(FV) i: " << i << " gradQ[i]: " << gradQ[i] << std::endl;
  }
  std::cout << std::endl;*/

  kernels::idx2 idx_gradQ(DIMENSIONS, NumberOfVariables+NumberOfParameters);

  //std::cout << idx_gradQ(0, 1) << " " << gradQ[idx_gradQ(0, 1)] << std::endl;
  //std::cout << idx_gradQ(1, 2) << " " << gradQ[idx_gradQ(1, 2)] << std::endl;
  //std::cout << std::endl;

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

  //const double mu = std::max(1.0e-10, Q[0]);
  //const double e = std::max(0.0, Q[3])/mu;
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> massGrad;
  massGrad(0) = gradQ[idx_gradQ(0, 0)];
  massGrad(1) = gradQ[idx_gradQ(1, 0)];

  Eigen::Matrix<double, MacroscaleInformation::dim, 1> energyGrad;
  energyGrad(0) = (gradQ[idx_gradQ(0, 3)] - energy*massGrad(0))/massDensity;
  energyGrad(1) = (gradQ[idx_gradQ(1, 3)] - energy*massGrad(1))/massDensity;

  Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> velJac;
  // momentum jacobian
  velJac(0, 0) = gradQ[idx_gradQ(0, 1)]; velJac(0, 1) = gradQ[idx_gradQ(1, 1)];
  velJac(1, 0) = gradQ[idx_gradQ(0, 2)]; velJac(1, 1) = gradQ[idx_gradQ(1, 2)];
  velJac = (velJac - velocity*massGrad.transpose())/massDensity;

  Eigen::Matrix<double, MacroscaleInformation::dim, 1> pressureGrad = energyGrad - velJac*velocity;
  const double pressureGradNorm = pressureGrad.norm();

  //std::cout << "pressure grad: " << pressureGrad.transpose() << std::endl;

  //const double mugradnorm = std::sqrt(gradQ[idx_gradQ(0, 0)]*gradQ[idx_gradQ(0, 0)]+gradQ[idx_gradQ(1, 0)]*gradQ[idx_gradQ(1, 0)]);
  const double kinetic = 0.5*velocity.dot(velocity);
  const double pressure = std::max(1.0e-10, energy-kinetic);

  const double divMomentum = gradQ[idx_gradQ(0, 1)] + gradQ[idx_gradQ(1, 2)];

  const double source = std::min(0.0, divMomentum);
  //const double source = divMomentum*std::abs(divMomentum);

  //std::cout << "energy: " << energy << " pressure: " << pressure << std::endl;

  //std::cout << "div momentum: " << divMomentum << std::endl;

  //std::cout << "scale: " << std::exp(-kinetic) << std::endl;

  //std::cout << "mass grad: " << massGrad.transpose() << std::endl;
  //std::cout << 1.0/(1.0+std::exp(-100.0*(massGrad.dot(massGrad)-1.0))) << std::endl;
  //std::cout << pressure*(1.0-1.0/(1.0+std::exp(-10.0*(massGrad.dot(massGrad)-10.0)))) << std::endl;
  //std::cout << std::endl;

  /*if( 1.0/(1.0+std::exp(-100.0*(pressureGrad.dot(pressureGrad)-5.0)))>0.5 ) {
    std::cout << "pressure grad: " << pressureGrad.transpose() << std::endl;
    std::cout << 1.0/(1.0+std::exp(-100.0*(pressureGrad.dot(pressureGrad)-5.0))) << std::endl;
  }*/


  BgradQ[0] = 0.0;
  BgradQ[1] = 0.0;
  BgradQ[2] = 0.0;
  BgradQ[3] = 0.0;
  //BgradQ[3] += 7.5*massDensity*energy*massGrad.transpose()*velJac*massGrad;
  //BgradQ[3] += 7.5*massDensity*(pressure/energy-0.9)*massGrad.dot(velocity);
  //BgradQ[3] += 10.0*massDensity/pressure*(velocity(0)*pressureGrad(0)+velocity(1)*pressureGrad(1));
  //BgradQ[3] += massDensity*energy*massGrad.dot(velocity);

  //BgradQ[3] = -kinetic/energy*massGrad.dot(massGrad);

  //BgradQ[3] = massDensity*kinetic*divMomentum;
  //BgradQ[3] = massGrad.dot(massGrad)*(velJac(0, 0)+velJac(1, 1))/std::max(1.0e-10, pressure);
  //BgradQ[3] = massDensity*( 1.0/(1.0+std::exp(-(std::pow(massGrad.norm()/std::max(1.0e-10, pressure), 2.0)-1.0))) );
  //BgradQ[3] = massGrad(0)-pressureGrad(0);
  //BgradQ[3] = massDensity*(pressureGrad(0) - pressureGrad(1)) + pressure*(massGrad(0) + massGrad(1));
  //BgradQ[3] = -pressure*massDensity*(massGrad(0) + massGrad(1));
  //BgradQ[3] = -3.0*pressure*massDensity*massDensity*source;
  //BgradQ[3] = pressure*massDensity*massDensity*(2.0/(1.0+std::exp(-1000.0*(pressureGrad.dot(pressureGrad)-0.5))));
  //BgradQ[3] = pressure*(2.0/(1.0+std::exp(-1000.0*(energyGrad.dot(energyGrad)-0.25))) + source);
  //BgradQ[3] = pressure*(1.0-1.0/(1+std::exp(-10.0*(massGrad.dot(massGrad)-10.0))));
  //BgradQ[3] = -massDensity*massDensity*pressure;
  //BgradQ[3] = e*mugradnorm*(1.0 + (gradQ[idx_gradQ(0, 1)] + gradQ[idx_gradQ(1, 2)] - Q[1]*gradQ[idx_gradQ(0, 0)]/Q[0] - Q[2]*gradQ[idx_gradQ(1, 0)]/Q[0]));
  //BgradQ[3] = 3.0*(gradQ[idx_gradQ(0, 0)]+gradQ[idx_gradQ(1, 0)]);

  //BgradQ[3] += energy*massDensity*divMomentum;
}
