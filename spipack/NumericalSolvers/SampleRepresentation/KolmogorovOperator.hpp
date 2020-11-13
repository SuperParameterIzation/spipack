#ifndef KOLMOGOROVOPERATOR_HPP_
#define KOLMOGOROVOPERATOR_HPP_

#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

namespace spi {
namespace NumericalSolvers {

/// Discretely represent the Kolmogorov operator \f$\mathcal{L}_{\psi,\alpha}\f$ using samples \f$\{\boldsymbol{x}^{(i)}\}\f$ from the distribution \f$\psi\f$
/**
Define the Kolmogorov operator applied to a smooth function \f$f\f$
\f{equation*}{
\mathcal{L}_{\psi,\alpha} f = \Delta f + (2-2\alpha) \nabla f \cdot \frac{\nabla \psi}{\psi}
\f}
Special cases:
- \f$\alpha=\frac{1}{2}\f$: The weighted Laplacian operator \f$\mathcal{L}_{\psi,0.5} f = \Delta f + \nabla f \cdot \frac{\nabla \psi}{\psi} = \psi^{-1} \nabla \cdot (\psi \nabla f) = \Delta_{\psi} f\f$
- \f$\alpha=1\f$: The Laplacian operator \f$\mathcal{L}_{\psi,1} f = \Delta f\f$

References:
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520315000020">"Variable bandwidth diffusion kernels" by T. Berry & J. Harlim</a>
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520317300982">"Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis</a>
*/
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

  /// Construct the Kolmogorov operator given samples from the underlying distribution \f$\psi\f$
  /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] options Setup options
  */
  KolmogorovOperator(std::shared_ptr<const spi::Tools::NearestNeighbors> const& samples, YAML::Node const& options);

  virtual ~KolmogorovOperator() = default;

  /// Tune the density estimate
  void TuneDensityEstimation();
private:

  /// Estimate the density of the underlying distribution
  DensityEstimation density;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
