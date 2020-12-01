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

In addition to the parameters/options below, this class has the same parameters/options as spi::NumericalSolvers::SampleRepresentation.

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"DensityOptions"   | <tt>YAML::Node</tt> | Same options as for the Kolmogorov operator | The options for the density estimation (see spi::NumericalSolvers::DensityEstimation) |
"OperatorParameter"   | <tt>double</tt> | <tt>1.0</tt> | The constant \f$c\f$ that defines the Kolmogorov operator |
"BandwidthExponent"   | <tt>double</tt> | <tt>-0.5</tt> | The numerical parameter \f$\beta\f$ used to discretize the Kolmogorov operator |
"NumEigenvalues"   | <tt>std::size_t</tt> | <tt>10</tt> | Store the smallest \f$n\f$ eigenvalues of the (symmetric) matrix \f$\boldsymbol{\hat{L}}\f$ (and their corresponding eigenvectors) |
"EigensolverTolerance"   | <tt>double</tt> | <tt>10^{-6}</tt> | The tolerance for the eigensolver |
"EigensolverMaxIterations"   | <tt>std::size_t</tt> | <tt>10^{3}</tt> | The maximum number of iterations for the eigensolver |

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
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* tune = &tuneDefault) override;

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
  @param[in] dens An estimate of the density function \f$\psi\f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}\f$
  \return The row sum of the kernel matrix
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::Ref<Eigen::MatrixXd> kmat) override;

  /// Construct the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  This function first estimates the density function \f$\psi_i\f$ at each sample \f$\boldsymbol{x}^{(i)}\f$.

  Given this density function, call another KernelMatrix function to fill the kernel matrix.

  @param[in] eps The bandwidth parameer \f$\epsilon\f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}\f$
  @param[in] tune <tt>true</tt>: Tune the parameter used to estimate the density \f$\psi\f$; <tt>false</tt>: Do not tune the parameter used to estimate \f$\psi\f$
  \return The row sum of the kernel matrix
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* tune = &tuneDefault) override;

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
  @param[in] dens An estimate of the density function \f$\psi\f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}\f$
  \return The row sum of the kernel matrix
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& dens, Eigen::SparseMatrix<double>& kmat) override;

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

  /// Get the number of eigenvalues/vectors of the (symmetric) matrix \f$\boldsymbol{\hat{L}}\f$ that we are storing
  /**
  \return The number of eigenvalues that we are storing
  */
  std::size_t NumEigenvalues() const;

  /// Get the tolerance for the eigensolver
  /**
  \return The eigensolver tolerance
  */
  double EigensolverTolerance() const;

  /// Get the maximum number of iterations for the eigensolver
  /**
  \return The eigensolver maximum number of iterations
  */
  double EigensolverMaxIterations() const;

  /// Compute the eigendecomposition of the matrix \f$\boldsymbol{\hat{L}}\f$
  /**
  The matrix \f$\boldsymbol{\hat{L}}\f$ is related to the discrete Kolmogorov operator \f$\boldsymbol{L}\f$  by a similarity transformation such that \f$\boldsymbol{S} \boldsymbol{L} \boldsymbol{S}^{-1} = \boldsymbol{\hat{L}}\f$, where \f$\boldsymbol{S}\f$ is a (full rank) diagonal matrix.
  @param[out] S The diagonal matrix \f$\boldsymbol{S}^{-1}\f$
  @param[out] Sinv The inverse diagonal matrix \f$\boldsymbol{S}^{-1}\f$
  @param[out] eigenvalues The eigenvalues of \f$\boldsymbol{\hat{L}}\f$
  @param[out] eigenvectors The eigenvectors of \f$\boldsymbol{\hat{L}}\f$
  */
  void ComputeEigendecomposition(Eigen::Ref<Eigen::VectorXd> S, Eigen::Ref<Eigen::VectorXd> Sinv, Eigen::Ref<Eigen::VectorXd> eigenvalues, Eigen::Ref<Eigen::MatrixXd> eigenvectors);

  /// Compute the coefficients for the expansion of a function \f$f\f$ using the eigenvectors as a basis
  /**
  We are given the eigendecomposition of the symmetric matrix \f$\boldsymbol{\hat{L}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T\f$ such that the discrete Kolmogorov operator is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{L}} \boldsymbol{S}\f$. The corresponding eigendecomposition of \f$\boldsymbol{L}\f$ is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T \boldsymbol{S}\f$

  Given a function \f$f\f$, this function computes the coefficients \f$\boldsymbol{c} = \boldsymbol{\hat{Q}}^T \boldsymbol{S} \boldsymbol{f}\f$, where \f$f^{(i)} = f(\boldsymbol{x}^{(i)})\f$.
  @param[in] S The diagonal matrix \f$\boldsymbol{S}\f$
  @param[in] eigenvectors The eigenvectors \f$\boldsymbol{\hat{Q}}\f$ of \f$\boldsymbol{\hat{L}}\f$.
  @param[in] f We will compute the coefficients such that \f$\boldsymbol{f} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{c}\f$ and \f$f(\boldsymbol{x}^{(i)}) \approx f^{(i)}\f$
  \return The coefficients \f$\boldsymbol{c}\f$
  */
  Eigen::VectorXd FunctionRepresentation(Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, std::function<double(Eigen::VectorXd const&)> f) const;

  /// Compute the coefficients for the expansion of a function \f$f\f$ using the eigenvectors as a basis
  /**
  We are given the eigendecomposition of the symmetric matrix \f$\boldsymbol{\hat{L}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T\f$ such that the discrete Kolmogorov operator is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{L}} \boldsymbol{S}\f$. The corresponding eigendecomposition of \f$\boldsymbol{L}\f$ is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T \boldsymbol{S}\f$

  Given a function \f$f\f$, this function computes the coefficients \f$\boldsymbol{c} = \boldsymbol{\hat{Q}}^T \boldsymbol{S} \boldsymbol{f}\f$, where \f$f^{(i)} = f(\boldsymbol{x}^{(i)})\f$.
  @param[in] S The diagonal matrix \f$\boldsymbol{S}\f$
  @param[in] eigenvectors The eigenvectors \f$\boldsymbol{\hat{Q}}\f$ of \f$\boldsymbol{\hat{L}}\f$.
  @param[in] feval We will compute the coefficients such that \f$\boldsymbol{f} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{c}\f$ and \f$f(\boldsymbol{x}^{(i)}) \approx f^{(i)}\f$
  \return The coefficients \f$\boldsymbol{c}\f$
  */
  Eigen::VectorXd FunctionRepresentation(Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, Eigen::Ref<const Eigen::VectorXd> const& feval) const;

  /// Compute the coefficients for the expansion of a function \f$f\f$ using the eigenvectors as a basis
  /**
  We are given the eigendecomposition of the discrete Kolmogorov operator \f$\boldsymbol{L} = \boldsymbol{Q}_l \boldsymbol{\Lambda} \boldsymbol{Q}_r\f$, where \f$\boldsymbol{Q}_l\f$ and \f$\boldsymbol{Q}_r\f$ are the left and right eigenvectors.

  Given a function \f$f\f$, this function computes the coefficients \f$\boldsymbol{c} = \boldsymbol{Q}_r \boldsymbol{f}\f$, where \f$f^{(i)} = f(\boldsymbol{x}^{(i)})\f$.
  @param[in] eigenvectors The right eigenvectors \f$\boldsymbol{Q}_r\f$ of \f$\boldsymbol{L}\f$.
  @param[in] f We will compute the coefficients such that \f$\boldsymbol{f} = \boldsymbol{Q}_l \boldsymbol{c}\f$ and \f$f(\boldsymbol{x}^{(i)}) \approx f^{(i)}\f$
  \return The coefficients \f$\boldsymbol{c}\f$
  */
  Eigen::VectorXd FunctionRepresentation(Eigen::Ref<const Eigen::MatrixXd> const& eigenvectorsRight, std::function<double(Eigen::VectorXd const&)> f) const;

  /// Compute the coefficients for the expansion of a function \f$f\f$ using the eigenvectors as a basis
  /**
  We are given the eigendecomposition of the discrete Kolmogorov operator \f$\boldsymbol{L} = \boldsymbol{Q}_l \boldsymbol{\Lambda} \boldsymbol{Q}_r\f$, where \f$\boldsymbol{Q}_l\f$ and \f$\boldsymbol{Q}_r\f$ are the left and right eigenvectors.

  Given a function \f$f\f$, this function computes the coefficients \f$\boldsymbol{c} = \boldsymbol{Q}_r \boldsymbol{f}\f$, where \f$f^{(i)} = f(\boldsymbol{x}^{(i)})\f$.
  @param[in] eigenvectors The right eigenvectors \f$\boldsymbol{Q}_r\f$ of \f$\boldsymbol{L}\f$.
  @param[in] feval We will compute the coefficients such that \f$\boldsymbol{f} = \boldsymbol{Q}_l \boldsymbol{c}\f$ and \f$f(\boldsymbol{x}^{(i)}) \approx f^{(i)}\f$
  \return The coefficients \f$\boldsymbol{c}\f$
  */
  Eigen::VectorXd FunctionRepresentation(Eigen::Ref<const Eigen::MatrixXd> const& eigenvectorsRight, Eigen::Ref<const Eigen::VectorXd> const& feval) const;

  /// Compute the pseudo-inverse of the Kolmogorov operator given its eigendecomposition
  /**
  We are given the eigendecomposition of the symmetric matrix \f$\boldsymbol{\hat{L}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T\f$ such that the discrete Kolmogorov operator is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{L}} \boldsymbol{S}\f$. The corresponding eigendecomposition of \f$\boldsymbol{L}\f$ is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T \boldsymbol{S}\f$.

  The smallest eigenvalue is \f$0\f$, so this function applies the pseudo-inverse to the right hand side function. We then translate the result so that its expected value is \f$0\f$.
  @param[in] rhs The right hand side function \f$\boldsymbol{r}\f$ evaluated at each sample
  @param[in] S The diagonal matrix \f$\boldsymbol{S}\f$
  @param[in] Sinv The inverse diagonal matrix \f$\boldsymbol{S}^{-1}\f$
  @param[in] eigenvalues The eigenvalues \f$\boldsymbol{\Lambda}\f$ OR the inverse eigenvalues \f$\boldsymbol{\Lambda}^{-1}\f$
  @param[in] eigenvectors The eigenvectors \f$\boldsymbol{\hat{Q}}\f$
  @param[in] inv <tt>true</tt>: The eigenvalue vector stores in the inverse eigenvalues, <tt>false</tt> (default): the eigenvalue vector stores the eigenvalues
  \return The function after applying the pseudo-inverse operator evaluated at each sample
  */
  Eigen::VectorXd PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& rhs, Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::VectorXd> const& Sinv, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, bool const inv = false) const;

  /// Compute the pseudo-inverse of the Kolmogorov operator given its eigendecomposition
  /**
  We are given the eigendecomposition of the symmetric matrix \f$\boldsymbol{\hat{L}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T\f$ such that the discrete Kolmogorov operator is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{L}} \boldsymbol{S}\f$. The corresponding eigendecomposition of \f$\boldsymbol{L}\f$ is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T \boldsymbol{S}\f$.

  The smallest eigenvalue is \f$0\f$, so this function applies the pseudo-inverse to the right hand side function. This function returns the coefficients of the expansion of the function after applying the pseudo-inverse. The function can be recovered as the weighed sum of the eigenvectors---it will not have expected value \f$0\f$.
  @param[in] rhs The right hand side function  \f$\boldsymbol{r}\f$ evaluated at each sample
  @param[in] S The diagonal matrix \f$\boldsymbol{S}\f$
  @param[in] eigenvalues The eigenvalues \f$\boldsymbol{\Lambda}\f$ OR the inverse eigenvalues \f$\boldsymbol{\Lambda}^{-1}\f$
  @param[in] eigenvectors The eigenvectors \f$\boldsymbol{\hat{Q}}\f$
  @param[in] inv <tt>true</tt>: The eigenvalue vector stores in the inverse eigenvalues, <tt>false</tt> (default): the eigenvalue vector stores the eigenvalues
  \return The function after applying the pseudo-inverse operator evaluated at each sample
  */
  Eigen::VectorXd PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& rhs, Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors, bool const inv = false) const;

  /// Compute the pseudo-inverse of the Kolmogorov operator given its eigendecomposition
  /**
  We are given the eigendecomposition of the symmetric matrix \f$\boldsymbol{\hat{L}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T\f$ such that the discrete Kolmogorov operator is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{L}} \boldsymbol{S}\f$. The corresponding eigendecomposition of \f$\boldsymbol{L}\f$ is \f$\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T \boldsymbol{S}\f$.

  The smallest eigenvalue is \f$0\f$, so this function applies the pseudo-inverse to the right hand side function. However, in this case we assume that instead of the right hand side function at each sample, we have the coeffients of this function for its expansion in the eigen basis. This function returns the coefficients of the expansion of the function after applying the pseudo-inverse. The function can be recovered as the weighed sum of the eigenvectors---it will not have expected value \f$0\f$.

  This function also returns the pseudo-inverse of the diagonal matrix with the eigenvalues as the diagonal. Storing this could save us the trouble of inverting the eigenvalues in the future.
  @param[in] rhs The coeffients of the right hand side function \f$\boldsymbol{r}\f$
  @param[in] S The diagonal matrix \f$\boldsymbol{S}\f$
  @param[in] eigenvalues The eigenvalues \f$\boldsymbol{\Lambda}\f$
  @param[in] eigenvectors The eigenvectors \f$\boldsymbol{\hat{Q}}\f$
  @param[out] eigenvaluesInv The pseudo-inverse eigenvalues \f$\boldsymbol{\Lambda}^{-1}\f$ (if the eigenvalues are zero, the pseudo-inverse eigenvalues are also zero)
  \return The coeffients of the function after applying the pseudo-inverse operator
  */
  Eigen::VectorXd PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& rhs, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, bool const inv = false) const;

  /// Compute the pseudo-inverse of the eigenvalues
  /**
  The pseudo inverse of the eigenvalues is the inverse, but the zero eigenvalue remains zero.
  @param[in] eigenvalues The eigenvalues
  \return The pseudo-inverse eigenvalues
  */
  Eigen::VectorXd PseudoInverse(Eigen::Ref<const Eigen::VectorXd> const& eigenvalues) const;

  /// Compute the gradient of a function represented in the eigenbasis
  /**
  */
  Eigen::MatrixXd FunctionGradient(Eigen::Ref<const Eigen::VectorXd> const& coeff, Eigen::Ref<const Eigen::VectorXd> const& S, Eigen::Ref<const Eigen::VectorXd> const& Sinv, Eigen::Ref<const Eigen::VectorXd> const& eigenvalues, Eigen::Ref<const Eigen::MatrixXd> const& eigenvectors) const;

private:

  /// The diagonal of the matrix \f$\boldsymbol{P}\f$, where \f$P^{(ii)} = \psi^{\beta}(\boldsymbol{x}^{(i)})\f$
  /**
  \f$\beta\f$ is the numerical parameter spi::NumericalSolvers::KolmogorovOperator::exponentPara.
  */
  Eigen::VectorXd P;

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

  /// The number of eigenvalues (and eigenvectors) of the (symmetric) matrix \f$\boldsymbol{\hat{L}}\f$ to store
  const std::size_t neig;

  /// The tolerance for the eigensolver
  const double eigensolverTol;

  /// The number of iterations for the eigensolver
  const std::size_t eigensolverMaxIt;

  /// The default values for the spi::NumericalSolvers::DensityEstimation class.
  struct DefaultParameters {
    /// The default operator parameter is \f$c=1\f$
    inline static const double operatorConstant = 1.0;

    /// The exponent parameter defaults to \f$\beta=-0.5\f$
    inline static const double exponentPara = -0.5;

    /// The default number of stored eigenvalues/vectors is \f$10\f$.
    inline static const std::size_t neig = 10;

    /// The default tolerance for the eigensolver is \f$10^{-6}\f$
    inline static const double eigensolverTol = 1.0e-6;

    /// The default number of iterations for the eigensolver is \f$10^{3}\f$
    inline static const std::size_t eigensolverMaxIt = 1e3;

  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
