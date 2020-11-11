#ifndef DENSITYESTIMATION_HPP_
#define DENSITYESTIMATION_HPP_

#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

namespace spi {
namespace NumericalSolvers {

/// Estimate the density of a probability distribution \f$\psi\f$ given samples \f$\{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$ such that \f$\boldsymbol{x}^{(i)} \sim \psi\f$.
/**
Let \f$\boldsymbol{K}_{\epsilon}\f$ be the kernel matrix computing in spi::NumericalSolvers::SampleRepresentation::KernelMatrix using the bandwidth parameter \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$.

In addition to the parameters/options below, this class has the same parameters/options as spi::NumericalSolvers::SampleRepresentation.

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"BandwidthParameter"   | <tt>double</tt> | <tt>1.0</tt> | The parmeter \f$\epsilon\f$ used to compute the kernel |
*/
class DensityEstimation : public SampleRepresentation {
public:

  /// Construct the density estimation by sampling a random variable from \f$\psi\f$
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] options Setup options
  */
  DensityEstimation(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node const& options);

  /// Construct the density estimation given samples from the underlying distribution \f$\psi\f$
  /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] options Setup options
  */
  DensityEstimation(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options);

  virtual ~DensityEstimation() = default;

  /// Get the bandwith parameter \f$\epsilon\f$
  /**
  \return The bandwith parameter \f$\epsilon\f$
  */
  double BandwidthParameter() const;

  /// Estimate the density at each sample
  /**
  \return The density estimation at each sample \f$\psi_i \approx \psi(\boldsymbol{x}^{(i)})\f$
  */
  Eigen::VectorXd Estimate() const;

private:

  /// The bandwidth parameter \f$\epsilon\f$
  const double bandwidthPara;

  /// The dimension of the manifold \f$m\f$
  const double manifoldDim;

  /// The default values for the spi::NumericalSolvers::DensityEstimation class.
  struct DefaultParameters {
    /// The default bandwidth parameter \f$\epsilon\f$ is \f$1.0\f$
    inline static const double bandwidthPara = 1.0;

    /// The default manifold dimension is \f$2\f$
    inline static const double manifoldDim = 2.0;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
