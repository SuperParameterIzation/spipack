#ifndef SAMPLEREPRESENTATION_HPP_
#define SAMPLEREPRESENTATION_HPP_

#include "spipack/Tools/NearestNeighbors.hpp"

namespace spi {
namespace NumericalSolvers {

/// Represent a distribution \f$\phi\f$ using samples
class SampleRepresentation {
public:

  /// Construct the sample representation by sampling a random variable from \f$\psi\f$
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] options Setup options
  */
  SampleRepresentation(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node const& options);

  virtual ~SampleRepresentation() = default;

private:
};

} // namespace NumericalSolvers
} // namespace spi

#endif
