#ifndef BOLTZMANN_HPP_
#define BOLTZMANN_HPP_

#include <sstream>

#include "spipack/SPIPACKDirectories.hpp"

#include "spipack/KineticEquations/KineticModels.hpp"

#include "spipack/ConservationEquations/ExaHyPE/Boltzmann/BoltzmannSolver.h"

// First, define your list
#define SPIPACK_BOLTZMANN_SOLVERMODES(F) \
  F(CompressibleEulerLimit), \
  F(ParticleMethods)

namespace spi {
namespace ConservationEquations {

class Boltzmann {
public:
  /// Get the exahype boltzmann solver
  static spiEX_Boltzmann::BoltzmannSolver* BoltzmannSolver();

  static void Initialize();

  /// Which mode are we using to solve the Boltzmann equation
  /**
  For example, we could solve the compressible Euler or represent the conditional velocity distribution using a particle method
  */
  #define F(e) e
  enum SolverMode { SPIPACK_BOLTZMANN_SOLVERMODES(F), NumSolverModes };
  #undef F

  /**
  @param[in] microscaleModels The micro-scale models that allow us to compute expectations with respect to the conditional velocity distribution \f$\nu(\boldsymbol{v} \vert \boldsymbol{x})\f$
  @param[in] options Numerical options for the Boltzmann model
  */
  Boltzmann(std::shared_ptr<spi::KineticEquations::KineticModels> const& microscaleModels, YAML::Node const& options);

  virtual ~Boltzmann() = default;

  /// Which mode are we running?
  /**
  \return The solver mode.
  */
  Boltzmann::SolverMode GetSolverMode() const;

  void Run();

private:
  /// A map from string to enum options for the solver mode
  static std::unordered_map<std::string, Boltzmann::SolverMode> SolverModesMap();

  /// Which mode are we running?
  const SolverMode mode;

  /// The file that defines the options for the exahype simulation
  const std::string constructFile = spi::SPIPACK_INSTALL_DIR+"/ExaHyPE/GenerateBoltzmann.exahype";
};

} // namespace ConservationEquations
} // namespace spi

#endif
