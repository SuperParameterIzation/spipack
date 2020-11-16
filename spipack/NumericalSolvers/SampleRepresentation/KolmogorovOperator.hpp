#ifndef KOLMOGOROVOPERATOR_HPP_
#define KOLMOGOROVOPERATOR_HPP_

#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

namespace spi {
namespace NumericalSolvers {

/// Discretely represent the Kolmogorov operator \f$\mathcal{L}_{\psi,c}\f$ using samples \f$\{\boldsymbol{x}^{(i)}\}\f$ from the distribution \f$\psi\f$
/**
Define the Kolmogorov operator applied to a smooth function \f$f\f$
\f{equation*}{
\mathcal{L}_{\psi,\alpha} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi}
\f}
Special cases:
- \f$c=0\f$: The Laplacian operator \f$\mathcal{L}_{\psi,0} f = \Delta f\f$
- \f$c=1\f$: The weighted Laplacian operator \f$\mathcal{L}_{\psi,1} f = \Delta f + \nabla f \cdot \frac{\nabla \psi}{\psi} = \psi^{-1} \nabla \cdot (\psi \nabla f) = \Delta_{\psi} f\f$

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

  /// Compute the density estimate
  /**
  @param[in] tune <tt>true</tt> (default): tune the bandwidth parameters, <tt>false</tt>: do not tune the bandwidth parameters (use stored values)
  \return The density estimation at each sample
  */
  Eigen::VectorXd EstimateDensity(bool const tune = true) const;

  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* tune = &tuneDefault) const override;

  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::Ref<Eigen::MatrixXd> kmat) const override;

  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* tune = &tuneDefault) const override;

  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::SparseMatrix<double>& kmat) const override;

private:

  /// By default, do we want to tune the bandwidth parameter values?
  /**
  The default value is <tt>true</tt>.
  */
  inline static const bool tuneDefault = true;

  /// Estimate the density of the underlying distribution
  std::shared_ptr<DensityEstimation> density;

  /// The operator parameter \f$\alpha\f$
  const double operatorConstant;

  /// Variable bandwidth exponent \f$\beta\f$
  const double variableBandwidthExponent;

  /// The default values for the spi::NumericalSolvers::DensityEstimation class.
  struct DefaultParameters {
    /// The default operator parameter is \f$c=1\f$
    inline static const double operatorConstant = 1.0;

    /// The default variable bandwidth exponent is \f$\beta=-0.5\f$
    inline static const double variableBandwidthExponent = -0.5;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
