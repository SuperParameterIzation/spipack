#ifndef KOLMOGOROVOPERATOR_HPP_
#define KOLMOGOROVOPERATOR_HPP_

#include "spipack/NumericalSolvers/SampleRepresentation/DensityEstimation.hpp"

namespace spi {
namespace NumericalSolvers {

/// Discretely represent the Kolmogorov operator \f$\mathcal{L}_{\psi,c}\f$ using samples \f$\{\boldsymbol{x}^{(i)}\}\f$ from the distribution \f$\psi\f$
/**
Define the Kolmogorov operator applied to a smooth function \f$f\f$
\f{equation*}{
\mathcal{L}_{\psi,c} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi}
\f}
Special cases:
- \f$c=0\f$: The Laplacian operator \f$\mathcal{L}_{\psi,0} f = \Delta f\f$
- \f$c=1\f$: The weighted Laplacian operator \f$\mathcal{L}_{\psi,1} f = \Delta f + \nabla f \cdot \frac{\nabla \psi}{\psi} = \psi^{-1} \nabla \cdot (\psi \nabla f) = \Delta_{\psi} f\f$

References:
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520315000020">"Variable bandwidth diffusion kernels" by T. Berry & J. Harlim</a>
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520317300982">"Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis</a>
*/
class KolmogorovOperator : public SampleRepresentation, public std::enable_shared_from_this<KolmogorovOperator> {
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

    /// Compute the squared bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  /**
  This function does NOT compute the kd trees. It assumes that spi::NumericalSolvers::SampleRepresentation::BuildKDTrees has already been called.

  \return The squared bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  */
  Eigen::VectorXd SquaredBandwidth() const;

  /// Compute the average of the pair-wise kernel derivative evaluations \f$\bar{k}^{\prime} = \frac{d}{d \epsilon} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] \f$
  /**
  We (recursively) compute the average
  \f{equation*}{
  \bar{k}^{\prime} = \frac{d}{d \epsilon} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] = -n^{-2} \sum_{i,j=0}^{n} \left. \frac{d k}{d \theta} \right|_{\theta = \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j}} \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon^2 r_i r_j}
  \f}
  @param[in] eps The bandwidth parameter \f$\epsilon\f$
  @param[in] rvec The vector \f$\boldsymbol{r}\f$
  \return The kernel derivative average \f$\bar{k}^{\prime}\f$.
  */
  virtual double KernelDerivativeAverage(double const eps, Eigen::VectorXd const& rvec) const override;

  /// Compute the average of the pair-wise kernel second derivative evaluations \f$\bar{k}^{\prime} = \frac{d^2}{d \epsilon^2} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] \f$
  /**
  We (recursively) compute the average
  \f{equation*}{
  \bar{k}^{\prime} = \frac{d^2}{d \epsilon^2} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] = -n^{-2} \sum_{i,j=0}^{n} \left( -2 \left. \frac{d k}{d \theta} \right|_{\theta = \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j}} \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon^3 r_i r_j} - \left. \frac{d^2 k}{d \theta^2} \right|_{\theta = \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j}} \left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon^2 r_i r_j} \right)^2 \right)
  \f}
  @param[in] eps The bandwidth parameter \f$\epsilon\f$
  @param[in] rvec The vector \f$\boldsymbol{r}\f$
  \return The kernel derivative average \f$\bar{k}^{\prime}\f$.
  */
  virtual double KernelSecondDerivativeAverage(double const eps, Eigen::VectorXd const& rvec) const override;

  virtual Eigen::MatrixXd KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* tune = &tuneDefault) const override;

  virtual Eigen::MatrixXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::Ref<Eigen::MatrixXd> kmat) const override;

  virtual Eigen::MatrixXd KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* tune = &tuneDefault) const override;

  virtual Eigen::MatrixXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::SparseMatrix<double>& kmat) const override;

  void TuneBandwidthParameter();

  /// Get the parameter used as the variable bandwidth exponent
  /**
  \return The variable bandwidth exponent \f$\alpha\f$.
  */
  double VariableBandwidthExponent() const;

  /// Get the exponent parameter \f$\beta\f$
  /**
  \return The exponent parameter \f$\beta\f$.
  */
  double ExponentParameter() const;

private:

  /// By default, do we want to tune the bandwidth parameter values?
  /**
  The default value is <tt>true</tt>.
  */
  inline static const bool tuneDefault = true;

  /// Estimate the density of the underlying distribution
  std::shared_ptr<DensityEstimation> density;

  /// The operator parameter \f$c\f$
  const double operatorConstant;

  /// The exponent parameter \f$\beta\f$
  const double exponentPara;

  /// The default values for the spi::NumericalSolvers::DensityEstimation class.
  struct DefaultParameters {
    /// The default operator parameter is \f$c=1\f$
    inline static const double operatorConstant = 1.0;

    /// The exponent parameter defaults to \f$\beta=-0.5\f$
    inline static const double exponentPara = -0.5;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
