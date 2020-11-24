#include "BoltzmannSolver_FV.h"

#include "BoltzmannSolver_FV_Variables.h"


tarch::logging::Log Boltzmann::BoltzmannSolver_FV::_log( "Boltzmann::BoltzmannSolver_FV" );

void Boltzmann::BoltzmannSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "Boltzmann::BoltzmannSolver_FV.h".

  // @todo Please implement/augment if required
}

void Boltzmann::BoltzmannSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  if (tarch::la::equals(t,0.0)) {
    Variables vars(Q);

    vars.massDensity() = 1.0+x[0];
    vars.momentum(1.0, 0.0);
    vars.coordinates(x[0], x[1]);
  }
}

void Boltzmann::BoltzmannSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "Boltzmann::BoltzmannSolver_FV.h".
  // Tip: See header file "Boltzmann::AbstractBoltzmannSolver_FV.h" for toolkit generated compile-time
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // @todo Please implement/augment if required
  lambda[0] = 1.0;
  lambda[1] = 1.0;
  lambda[2] = 1.0;
  lambda[3] = 1.0;
  lambda[4] = 1.0;
}

void Boltzmann::BoltzmannSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  // Tip: You find documentation for this method in header file "Boltzmann::BoltzmannSolver_FV.h".
  // Tip: See header file "Boltzmann::AbstractBoltzmannSolver_FV.h" for toolkit generated compile-time
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // @todo Please implement/augment if required
  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];
  stateOutside[4] = stateInside[4];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void Boltzmann::BoltzmannSolver_FV::flux(const double* const Q,double** const F) {
  // Tip: You find documentation for this method in header file "Boltzmann::BoltzmannSolver_FV.h".
  // Tip: See header file "Boltzmann::AbstractBoltzmannSolver_FV.h" for toolkit generated compile-time
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.

  // @todo Please implement/augment if required
  F[0][0] = 0.0;
  F[0][1] = 0.0;
  F[0][2] = 0.0;
  F[0][3] = 0.0;
  F[0][4] = 0.0;

  F[1][0] = 0.0;
  F[1][1] = 0.0;
  F[1][2] = 0.0;
  F[1][3] = 0.0;
  F[1][4] = 0.0;

}




//You can either implement this method or modify fusedSource
void Boltzmann::BoltzmannSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  // Tip: You find documentation for this method in header file "Boltzmann::BoltzmannSolver_FV.h".
  // Tip: See header file "Boltzmann::AbstractBoltzmannSolver_FV.h" for toolkit generated compile-time
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  // @todo Please implement/augment if required
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}
