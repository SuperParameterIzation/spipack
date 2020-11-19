// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
// ==============================================
// Please do not change the implementations below
// ==============================================
#include "CompressibleEulerSolver.h"

#include "kernels/limiter/generic/Limiter.h"
#include "kernels/LimiterProjectionMatrices.h"

CompressibleEuler::CompressibleEulerSolver::CompressibleEulerSolver(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int haloBufferCells,
        const int limiterBufferCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int DMPObservables,
        const double DMPRelaxationParameter,
        const double DMPDifferenceScaling
) :
  exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
      "CompressibleEulerSolver",
    new CompressibleEuler::CompressibleEulerSolver_ADERDG(
      maximumMeshSize,maximumMeshDepth,haloCells,haloBufferCells,limiterBufferCells,regularisedFineGridLevels,timeStepping,DMPObservables),
    new CompressibleEuler::CompressibleEulerSolver_FV(
      maximumMeshSize, timeStepping),
    DMPRelaxationParameter,
    DMPDifferenceScaling) {
  kernels::computeDG2FVProjector<Order,PatchSize>(dg2fv);    
  kernels::computeFV2DGProjector<Order,PatchSize>(fv2dg,dg2fv);
  kernels::computeLegendre2LobattoProjector<Order>(leg2lob);  
}

void CompressibleEuler::CompressibleEulerSolver::projectOnFVLimiterSpace(const double* const luh, double* const lim) const {
  kernels::limiter::generic::c::projectOnFVLimiterSpace<NumberOfVariables+NumberOfParameters, Order+1, PatchSize, GhostLayerWidth>(luh, dg2fv, lim);
}

void CompressibleEuler::CompressibleEulerSolver::projectOnDGSpace(const double* const lim, double* const luh) const {
  kernels::limiter::generic::c::projectOnDGSpace<NumberOfVariables+NumberOfParameters, Order+1, PatchSize, GhostLayerWidth>(lim, fv2dg, luh);
}

bool CompressibleEuler::CompressibleEulerSolver::discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* const boundaryMinPerVariables, double* const boundaryMaxPerVariables) {
  return kernels::limiter::generic::c::discreteMaximumPrincipleAndMinAndMaxSearch<AbstractCompressibleEulerSolver_ADERDG, NumberOfDMPObservables, Order+1, PatchSize, GhostLayerWidth>(
    luh, *static_cast<AbstractCompressibleEulerSolver_ADERDG*>(_solver.get()), dg2fv, fv2dg, leg2lob, _DMPMaximumRelaxationParameter, _DMPDifferenceScaling, boundaryMinPerVariables, boundaryMaxPerVariables);
}

void CompressibleEuler::CompressibleEulerSolver::findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) {
  kernels::limiter::generic::c::findCellLocalMinAndMax<AbstractCompressibleEulerSolver_ADERDG, NumberOfDMPObservables, Order+1, PatchSize>(
    luh, *static_cast<AbstractCompressibleEulerSolver_ADERDG*>(_solver.get()), dg2fv, fv2dg, leg2lob, localMinPerVariables, localMaxPerVariable);
}

void CompressibleEuler::CompressibleEulerSolver::findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) {
  kernels::limiter::generic::c::findCellLocalLimiterMinAndMax<AbstractCompressibleEulerSolver_ADERDG, NumberOfDMPObservables, Order+1, PatchSize, GhostLayerWidth>(
    lim, *static_cast<AbstractCompressibleEulerSolver_ADERDG*>(_solver.get()), localMinPerObservable, localMaxPerObservable);
}