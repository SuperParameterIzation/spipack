#ifndef SAMPLEREPRESENTATION_HPP_
#define SAMPLEREPRESENTATION_HPP_

#include "spipack/Tools/NearestNeighbors.hpp"

namespace spi {
namespace NumericalSolvers {

/// Represent a distribution \f$\psi\f$ using samples
class SampleRepresentation {
public:

  /// Construct the sample representation by sampling a random variable from \f$\psi\f$
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] options Setup options
  */
  SampleRepresentation(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node const& options);

  /// Construct the sample representation given samples from the underlying distribution \f$\psi\f$
  /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] options Setup options
  */
  SampleRepresentation(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options);

  virtual ~SampleRepresentation() = default;

  /// How many samples are in the collection?
  /**
    \return The number of samples in the sample collection (GraphLaplacian::samples)
  */
  std::size_t NumSamples() const;

  /// Get the \f$i^{th}\f$ point from the point cloud.
  /**
  @param[in] i We want to get this point
  \return The f$i^{th}\f$ point from the point cloud.
  */
  Eigen::Ref<Eigen::VectorXd const> Point(std::size_t const i) const;

protected:

  /// Store the samples from \f$\psi\f$.
  const spi::Tools::NearestNeighbors samples;

private:
};

} // namespace NumericalSolvers
} // namespace spi

#endif
