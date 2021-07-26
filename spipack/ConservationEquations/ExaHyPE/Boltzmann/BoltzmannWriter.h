// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef POSTPROCESSING_BoltzmannWriter_CLASS_HEADER_
#define POSTPROCESSING_BoltzmannWriter_CLASS_HEADER_

#include "exahype/plotters/Plotter.h"

namespace spiEX_Boltzmann {
  class BoltzmannSolver;
  class BoltzmannWriter;
}

class spiEX_Boltzmann::BoltzmannWriter : public exahype::plotters::Plotter::UserOnTheFlyPostProcessing {
public:
  BoltzmannWriter(spiEX_Boltzmann::BoltzmannSolver& solver);
  virtual ~BoltzmannWriter();

  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp) override;
};

#endif /* POSTPROCESSING_BoltzmannWriter_CLASS_HEADER_ */