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
  Compute an estimate of the underlying \f$\psi\f$ at each sample using spi::NumericalSolvers::DensityEstimation.
  @param[in] tune <tt>true</tt> (default): tune the bandwidth parameters for the density estimation, <tt>false</tt>: do not tune the bandwidth parameters (use stored values)
  \return The density estimation at each sample
  */
  Eigen::VectorXd EstimateDensity(bool const tune = true) const;

  /// Construct the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  This function first estimates the density function \f$\psi_i\f$ at each sample \f$\boldsymbol{x}^{(i)}\f$.

  Given this density function, call another KernelMatrix function to fill the kernel matrix.

  @param[in] eps The bandwidth parameer \f$\epsilon\f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}\f$
  @param[in] tune <tt>true</tt>: Tune the parameter used to estimate the density \f$\psi\f$; <tt>false</tt>: Do not tune the parameter used to estimate \f$\psi\f$
  \return The row sum of the kernel matrix
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* tune = &tuneDefault) const override;

  /// Construct the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  This function first computes the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is
  \f{equation*}{
    \widetilde{K}_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon \psi_i^{\beta} \psi_j^{\beta}} \right).
  \f}
  It then computes the normalization
  \f{equation*}{
  \tilde{\rho}^{(i)} = \frac{1}{\psi_i^{\beta \alpha d}} \sum_{j=1}^{n} \widetilde{K}_{\epsilon}^{(ij)},
  \f}
  where \f$\alpha = 1 + \frac{1}{2} d \beta + \beta - \frac{1}{2} c\f$. We finally fill the kernel matrix with
  \f{equation*}{
    K_{\epsilon}^{(ij)} = \frac{\widetilde{K}_{\epsilon}^{(ij)}}{\tilde{\rho}^{(i)} \tilde{\rho}^{(j)}}.
  \f}
  and return the row sum vector
  \f{equation*}{
  \rho^{(i)} = \sum_{j=1}^{n} K_{\epsilon}^{(ij)},
  \f}
  @param[in] eps The bandwidth parameer \f$\epsilon\f$
  @param[in] dens An estimate of the density function \f$\psi|f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}\f$
  \return The row sum of the kernel matrix
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::Ref<Eigen::MatrixXd> kmat) const override;

  /// Construct the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  This function first estimates the density function \f$\psi_i\f$ at each sample \f$\boldsymbol{x}^{(i)}\f$.

  Given this density function, call another KernelMatrix function to fill the kernel matrix.

  @param[in] eps The bandwidth parameer \f$\epsilon\f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}\f$
  @param[in] tune <tt>true</tt>: Tune the parameter used to estimate the density \f$\psi\f$; <tt>false</tt>: Do not tune the parameter used to estimate \f$\psi\f$
  \return The row sum of the kernel matrix
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* tune = &tuneDefault) const override;

  /// Construct the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  This function first computes the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is
  \f{equation*}{
    \widetilde{K}_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon \psi_i^{\beta} \psi_j^{\beta}} \right).
  \f}
  It then computes the normalization
  \f{equation*}{
  \tilde{\rho}^{(i)} = \frac{1}{\psi_i^{\beta \alpha d}} \sum_{j=1}^{n} \widetilde{K}_{\epsilon}^{(ij)},
  \f}
  where \f$\alpha = 1 + \frac{1}{2} d \beta + \beta - \frac{1}{2} c\f$. We finally fill the kernel matrix with
  \f{equation*}{
    K_{\epsilon}^{(ij)} = \frac{\widetilde{K}_{\epsilon}^{(ij)}}{\tilde{\rho}^{(i)} \tilde{\rho}^{(j)}}.
  \f}
  and return the row sum vector
  \f{equation*}{
  \rho^{(i)} = \sum_{j=1}^{n} K_{\epsilon}^{(ij)},
  \f}
  @param[in] eps The bandwidth parameer \f$\epsilon\f$
  @param[in] dens An estimate of the density function \f$\psi|f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}\f$
  \return The row sum of the kernel matrix
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::SparseMatrix<double>& kmat) const override;

  /// Tune the bandwidth parameter for the Kolmogorov operator
  /**
  @param tuneDens Also tune the bandwidth for the underlying density estimation (defaults to <tt>true</tt>)
  */
  void TuneBandwidthParameter(bool const tuneDens = true);

  /// Get the parameter used as the variable bandwidth exponent \f$\alpha\f$
  /**
  \return The variable bandwidth exponent \f$\alpha = 1 + \frac{1}{2} d \beta + \beta - \frac{1}{2} c\f$.
  */
  double VariableBandwidthExponent() const;

  /// Get the exponent parameter \f$\beta\f$
  /**
  \return The exponent parameter \f$\beta\f$.
  */
  double ExponentParameter() const;

  /// Get the density estimator
  /**
  \return The density estimator
  */
  std::shared_ptr<DensityEstimation> Density() const;

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
