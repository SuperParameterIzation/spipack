#ifndef BOLTZMANN_HPP_
#define BOLTZMANN_HPP_

#include "exahype/runners/Runner.h"

namespace exahype {
  namespace parser {
    /// Forward declaration of exahype::parser::Parser class
    class Parser;
  }
}

namespace spi {
namespace ConservationEquations {

class Boltzmann : public exahype::runners::Runner {
public:
  Boltzmann();

  virtual ~Boltzmann() = default;

  void Run();

private:

  inline static exahype::parser::Parser parser;
  inline static std::vector<std::string> cmdlineargs;

  void Run(exahype::repositories::Repository& repository);

  /// Construct the <tt>ExaHyPE</tt> (external library) parser based on the internal configuration file
  /**
  \return parser The <tt>ExaHyPE</tt> parser
  */
  void ConstructBoltzmannParser(std::stringstream& specfile);

  void UpdateSmallScale(); 
};

} // namespace ConservationEquations
} // namespace spi

#endif
