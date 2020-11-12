#ifndef KOLMOGOROVOPERATOR_HPP_
#define KOLMOGOROVOPERATOR_HPP_

#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

namespace spi {
namespace NumericalSolvers {

class KolmogorovOperator : public SampleRepresentation {
public:

  /// Construct the Kolmogorov operator by sampling a random variable from \f$\psi\f$
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] options Setup options
  */
  KolmogorovOperator(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node const& options);

  /// Construct the Kolmogorov operator given samples from the underlying distribution \f$\psi\f$
  /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] options Setup options
  */
  KolmogorovOperator(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options);

  virtual ~KolmogorovOperator() = default;
private:
};

} // namespace NumericalSolvers
} // namespace spi

#endif
