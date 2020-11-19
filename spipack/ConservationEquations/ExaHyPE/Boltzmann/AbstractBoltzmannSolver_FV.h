// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#ifndef __AbstractBoltzmannSolver_FV_CLASS_HEADER__
#define __AbstractBoltzmannSolver_FV_CLASS_HEADER__

#include <ostream>

#include "exahype/solvers/FiniteVolumesSolver.h"

/**
 * We include Peano's assertion collection here.
 */
#include "tarch/Assertions.h"

namespace Boltzmann{
  class AbstractBoltzmannSolver_FV;
  class BoltzmannSolver_FV;
}

class Boltzmann::AbstractBoltzmannSolver_FV : public exahype::solvers::FiniteVolumesSolver {
  public:
    static constexpr int NumberOfVariables         = 5;
    static constexpr int NumberOfParameters        = 0;
    static constexpr int NumberOfGlobalObservables = 0;
    static constexpr int PatchSize                 = 7;
    static constexpr int GhostLayerWidth           = 2;
    static constexpr double CFL                    = 0.9;
    
    // virtual getters for the constexpr's ; TODO(Dominic): Who uses these methods?
    int constexpr_getNumberOfVariables()  const { return NumberOfVariables; };
    int constexpr_getNumberOfParameters() const { return NumberOfParameters; };
    int constexpr_getPatchSize()          const { return PatchSize; };
    int constexpr_getGhostLayerWidth()    const { return GhostLayerWidth; }
    double constexpr_getCFLNumber()       const { return CFL; };
  
    // Access Objects
    class VariableMetrics;
    class Variables;
    class ReadOnlyVariables;
    class Fluxes;
    class VariableShortcuts;
    class VariableMultiplicities;
    class VariableNames;
    class GlobalObservables;
    class ReadOnlyGlobalObservables;
    
    AbstractBoltzmannSolver_FV(
        const double maximumMeshSize,
        const exahype::solvers::Solver::TimeStepping timeStepping
);
    
    void solutionUpdate(double* luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize, const double t, const double dt,double& maxAdmissibleDt) override;
    
    double stableTimeStepSize(const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& cellSize) override;
    void adjustSolution(double* const luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,double t,double dt) override;
    
    void ghostLayerFilling(double* const luh,const double* const luhNeighbour,const tarch::la::Vector<DIMENSIONS,int>& neighbourPosition) override;
    void ghostLayerFillingAtBoundary(double* const luh,const double* const luhbnd,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) override;
    void boundaryLayerExtraction(double* const luhbnd,const double* const luh,const tarch::la::Vector<DIMENSIONS,int>& boundaryPosition) override;
    void boundaryConditions(double* const luh,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const tarch::la::Vector<DIMENSIONS, int>& posCell,const tarch::la::Vector<DIMENSIONS, int>& posBoundary) override;
    void resetGlobalObservables(double* const globalObservables) final override;
    void updateGlobalObservables(
        double* const                               globalObservables,
        const double* const                         luh,
        const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS,double>& cellSize,
        const double t,
        const double dt) final override;
    void mergeGlobalObservables(double* const globalObservables,const double* const otherObservables) final override;
    void wrapUpGlobalObservables(double* const globalObservables) final override;
    
    /// Apr 18, Coding Week: Riemann Solvers in FV. Hopefully inlined as evaluated point wise.
    double riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int direction) override;

    static void constantsToString(std::ostream& os);
    static void abortWithMsg(const char* const msg);
    
    /**
     * @defgroup User PDE kernels
     */
    ///@{
    /**
     * Compute the eigenvalues of the flux tensor per coordinate direction \p d.
     *
     * \param[in] Q  the conserved variables associated with a quadrature node
     *               as C array (already allocated).
     * \param[in] d  the column of the flux vector (d=0,1,...,DIMENSIONS).
     * \param[inout] lambda the eigenvalues as C array (already allocated).
     */
    virtual void eigenvalues(const double* const Q,const int direction,double* const lambda) {}
  
    /**
     * Impose boundary conditions at a point on a boundary face
     * within the time interval [t,t+dt].
     *
     * \param[in]    x         the physical coordinate on the face.
     * \param[in]    t         the start of the time interval.
     * \param[in]    dt        the width of the time interval.
     * \param[in]    faceIndex indexing of the face (0 -- {x[0]=xmin}, 1 -- {x[1]=xmax}, 2 -- {x[1]=ymin}, 3 -- {x[2]=ymax}, and so on,
     *                         where xmin,xmax,ymin,ymax are the bounds of the cell containing point x.
     * \param[in]    d         the coordinate direction the face normal is pointing to.
     * \param[in]    QIn       the conserved variables at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[in]    FIn       the normal fluxes at point x from inside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[inout] QOut      the conserved variables at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     * \param[inout] FOut      the normal fluxes at point x from outside of the domain
     *                         and time-averaged (over [t,t+dt]) as C array (already allocated).
     */
    virtual void boundaryValues(const double* const x,const double t,const double dt,const int faceIndex,const int direction,const double* const stateIn,double* const stateOut) {}
  
    /**
     * Compute a pointSource contribution.
     * 
     * @TODO: Document me, please.
     **/
     virtual void pointSource(const double* const Q,const double* const x,const double t,const double dt, double* const forceVector,int n) {}
  
    /**
     * Compute the Algebraic Sourceterms.
     * 
     * You may want to overwrite this with your PDE Source (algebraic RHS contributions).
     * However, in all schemes we have so far, the source-type contributions are
     * collected with non-conservative contributions into a fusedSource, see the
     * fusedSource method. From the kernels given with ExaHyPE, only the fusedSource
     * is called and there is a default implementation for the fusedSource calling
     * again seperately the nonConservativeProduct function and the algebraicSource
     * function.
     *
     * \param[in]    Q the conserved variables (and parameters) associated with a quadrature point
     *                 as C array (already allocated).
     * \param[inout] S the source point as C array (already allocated).
     */
    virtual void algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {}

    /**
     * Compute the nonconservative term $B(Q) \nabla Q$.
     * 
     * This function shall return a vector BgradQ which holds the result
     * of the full term. To do so, it gets the vector Q and the matrix
     * gradQ which holds the derivative of Q in each spatial direction.
     * Currently, the gradQ is a continous storage and users can use the
     * kernels::icellSize2 class in order to compute the positions inside gradQ.
     *
     * @TODO: Check if the following is still right:
     * 
     * !!! Warning: BgradQ is a vector of size NumberOfVariables if you
     * use the ADER-DG kernels for nonlinear PDEs. If you use
     * the kernels for linear PDEs, it is a tensor with dimensions
     * Dim x NumberOfVariables.
     * 
     * \param[in]   Q   the vector of unknowns at the given position
     * \param[in]   gradQ   the gradients of the vector of unknowns,
     *                  stored in a linearized array.
     * \param[inout]  The vector BgradQ (extends nVar), already allocated. 
     *
     **/
    virtual void nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {}
    
    /**
     * Compute the conserved flux.
     * 
     * \param[in]  Q the conserved variabels (and parameters) associated with a
     *               quadrature point as C array.
     * \param[inout] F a C array with shape [nDim][nVars]. That is, this is an C list
     *               holding pointers to actual lists. Thus, the storage may be noncontinous.
     *               In any case, the storage has already been allocated.
     **/
    virtual void flux(const double* const Q,double** const F) {}
  
    /**
     * Compute the flux with diffusive component.
     * 
     * \param[in]  Q the conserved variabels (and parameters) associated with a
     *               quadrature point as C array.
     *             gradQ the gradient of the conserved variabels (and parameters) associated
     *                   with a quadrature point as C array.
     * \param[inout] F a C array with shape [nDim][nVars]. That is, this is an C list
     *               holding pointers to actual lists. Thus, the storage may be noncontinous.
     *               In any case, the storage has already been allocated.
     **/
    virtual void viscousFlux(const double* const Q, const double* const gradQ, double** const F) {}
    
    virtual void viscousEigenvalues(const double* const Q,const int direction,double* const lambda) {}

    /**
     * Resets the vector of global observables to some suitable initial value, e.g.
     * the smallest possible double if one wants to compute the maximum.
     *
     *\param[inout] globalObservables The global observables that we want to reset.
     */
    virtual void resetGlobalObservables(GlobalObservables& globalObservables) const;
    
    /**
     * Computes observables from a cell's solution values.
     * The result will be merged with the global observables.
     *
     *\param[inout] globalObservables The mapped observables.
     *\param[in]    luh               The solution array.
     *\param[in]    cellSize          The size of a cell.
     */
    virtual void mapGlobalObservables(
        GlobalObservables&                          globalObservables,
        const double* const                         luh,
        const tarch::la::Vector<DIMENSIONS,double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS,double>& cellSize,
        const double t,
        const double dt) const;

    /**
     * This method merges two vectors of (global) observables.
     *
     * For example, if one wants to compute the maximum of global variable i,
     * one should set
     *
     * observables[i] = std::max( globalObservables[i], otherObservables[i])
     *
     * and so on.
     *
     *\param[inout] globalObservables  The (merged) observables.
     *\param[in]    otherObservables   other observables that we want to merge with the first argument.
     */
    virtual void mergeGlobalObservables(
        GlobalObservables&         globalObservables,
        ReadOnlyGlobalObservables& otherObservables) const;
    
    /**
     * Wrap up the global observables, e.g. finalise an L^2 integral
     * by applying a square root.
     *
     * @note This routine is only called on the global master, after
     * the reduction has completed.
     * 
     * @note If this FV solver is managed by a LimitingADERDGSolver,
     * only the wrap up method of the ADER-DG solver must be implemented
     * as the wrapped up observables will be copied from the ADER-DG to the FV solver. 
     *
     *\param[inout] observables  The global observables.
     */
    virtual void wrapUpGlobalObservables(GlobalObservables& globalObservables) const;
    
    /**
     * Signals a user solver that ExaHyPE just started a new time step.
     *
     * @param[in] minTimeStamp the minimum time stamp (over all cells)
     *
     * @note This function is invoked before the first predictor computation
     * when the non fused time stepping is run. Otherwise, it is invoked after
     * the first predictor computation. It will always be called before
     * "adjustSolution" is invoked.
     *
     * @note [MPI] For convenient reductions consider to specify "global_observables" in
     * the specification file.
     * 
     * @param isFirstTimeStepOfBatchOrNoBatch  if this is the first time step of a batch or no batch is run.
     */
    virtual void beginTimeStep(const double minTimeStamp,const bool isFirstTimeStepOfBatchOrNoBatch) {}

    /**
     * Signals a user solver that ExaHyPE just finished a time step.
     *
     * @param[in] minTimeStamp the minimum time stamp (over all cells)
     *
     * @note This function is invoked after the solution was updated in all
     * cells.
     *
     * @note [MPI] For convenient reductions consider to specify "global_observables" in
     * the specification file. 
     *
     * @param isFirstTimeStepOfBatchOrNoBatch  if this is the first time step of a batch or no batch is run.
     * @param isLastTimeStepOfBatchOrNoBatch   if this is the last time step of a batch or if no batch is run.
     */
    virtual void endTimeStep(const double minTimeStamp,const bool isLastTimeStepOfBatchOrNoBatch)   {}
};


#endif // __AbstractBoltzmannSolver_FV_CLASS_HEADER__