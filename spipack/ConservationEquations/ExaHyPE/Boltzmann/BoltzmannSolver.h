#ifndef __BoltzmannSolver_CLASS_HEADER__
#define __BoltzmannSolver_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
//
// ========================
//   www.exahype.eu
// ========================

#include <string>

#include "exahype/solvers/LimitingADERDGSolver.h"
#include "BoltzmannSolver_ADERDG.h"
#include "BoltzmannSolver_FV.h"

#include "spipack/KineticEquations/KineticModels.hpp"

namespace spiEX_Boltzmann{
  class BoltzmannSolver;
}

class spiEX_Boltzmann::BoltzmannSolver: public exahype::solvers::LimitingADERDGSolver {
  public:
    static constexpr int NumberOfVariables         = spiEX_Boltzmann::AbstractBoltzmannSolver_ADERDG::NumberOfVariables;
    static constexpr int NumberOfParameters        = spiEX_Boltzmann::AbstractBoltzmannSolver_ADERDG::NumberOfParameters;
    static constexpr int Order                     = spiEX_Boltzmann::AbstractBoltzmannSolver_ADERDG::Order;
    static constexpr int NumberOfGlobalObservables = spiEX_Boltzmann::AbstractBoltzmannSolver_ADERDG::NumberOfGlobalObservables;
    static constexpr int NumberOfDMPObservables    = spiEX_Boltzmann::AbstractBoltzmannSolver_ADERDG::NumberOfDMPObservables;
    static constexpr int PatchSize                 = spiEX_Boltzmann::AbstractBoltzmannSolver_FV::PatchSize;
    static constexpr int GhostLayerWidth           = spiEX_Boltzmann::AbstractBoltzmannSolver_FV::GhostLayerWidth;

    // limiter projection matrices
    double dg2fv[(Order+1)*PatchSize];
    double fv2dg[(Order+1)*PatchSize];
    double leg2lob[(Order+1)*(Order+1)];

    BoltzmannSolver(
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
);

    void projectOnFVLimiterSpace(const double* const luh, double* const lim) const override;
    void projectOnDGSpace(const double* const lim, double* const luh) const override;
    bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* const boundaryMinPerVariables, double* const boundaryMaxPerVariables) override;
    void findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) override;
    void findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) override;

    /// The kinetic models, used to update the micro-scale simulations
    std::shared_ptr<spi::KineticEquations::KineticModels> kineticModels;
};

#endif // __BoltzmannSolver_CLASS_HEADER__
