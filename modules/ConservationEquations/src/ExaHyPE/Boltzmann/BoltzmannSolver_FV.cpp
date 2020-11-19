#include "BoltzmannSolver_FV.h"

#include "BoltzmannSolver_FV_Variables.h"


tarch::logging::Log Boltzmann::BoltzmannSolver_FV::_log( "Boltzmann::BoltzmannSolver_FV" );

void Boltzmann::BoltzmannSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
  // Tip: You find documentation for this method in header file "Boltzmann::BoltzmannSolver_FV.h".
  
  // @todo Please implement/augment if required
}

void Boltzmann::BoltzmannSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
  // Tip: You find documentation for this method in header file "Boltzmann::BoltzmannSolver_FV.h".
  // Tip: See header file "Boltzmann::AbstractBoltzmannSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  Q[0] = 0.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = 0.0;
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




void Boltzmann::BoltzmannSolver_FV::viscousFlux(const double* const Q,const double* const gradQ, double** const F) {
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

void Boltzmann::BoltzmannSolver_FV::viscousEigenvalues(const double* const Q, const int dIndex, double* const lambda) {
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




