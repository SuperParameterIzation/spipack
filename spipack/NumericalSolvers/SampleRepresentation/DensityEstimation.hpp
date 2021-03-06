#ifndef DENSITYESTIMATION_HPP_
#define DENSITYESTIMATION_HPP_

#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

namespace spi {
namespace NumericalSolvers {

/// Estimate the density of a probability distribution \f$\psi\f$ given samples \f$\{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$ such that \f$\boldsymbol{x}^{(i)} \sim \psi\f$.
/**
Let \f$\boldsymbol{K}_{\epsilon}\f$ be the kernel matrix computing in spi::NumericalSolvers::SampleRepresentation::KernelMatrix using the bandwidth parameter \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$. We estimate the density at each sample as
\f{equation*}{
    \psi^{(i)}_{\epsilon} = \sum_{j=1}^{n} \frac{K_{\epsilon}^{(ij)}}{n (\pi \epsilon r_i^2)^{m/2}},
\f}
where \f$m\f$ is the manifold dimension.

In addition to the parameters/options below, this class has the same parameters/options as spi::NumericalSolvers::SampleRepresentation.

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"TuneManifoldDimension"   | <tt>bool</tt> | <tt>false</tt> | Tune the manifold dimension \f$m\f$; if we know the exact manifold dimension then we do not need to tune it. |

References:
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520315000020">"Variable bandwidth diffusion kernels" by T. Berry & J. Harlim</a>
- <a href="https://www.sciencedirect.com/science/article/pii/S1063520317300982">"Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis</a>
*/
class DensityEstimation : public SampleRepresentation, public std::enable_shared_from_this<DensityEstimation> {
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

  /// Construct the density estimation given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  DensityEstimation(std::shared_ptr<const spi::Tools::NearestNeighbors> const& samples, YAML::Node const& options);

  virtual ~DensityEstimation() = default;

  /// Estimate the density at each sample
  /**
  This function computes the bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$ and the kernel marix \f$\boldsymbol{K}_{\epsilon}\f$. We then estimate the density at each sample as
  \f{equation*}{
  \psi^{(i)}_{\epsilon} = \sum_{j=1}^{n} \frac{K_{\epsilon}^{(ij)}}{n (\pi \epsilon r_i^2)^{m/2}},
  \f}
  where \f$m\f$ is the manifold dimension.

  This function does NOT compute the kd trees. It assumes that spi::NumericalSolvers::SampleRepresentation::BuildKDTrees has already been called.

  If <tt>tune</tt> is <tt>true</tt>, we must choose the optimal \f$\epsilon\f$ for the density estimation. Given the bandwidth exponent \f$l \in [l_{min}, l_{max}]\f$, we choose \f$\epsilon \in [2^{(l_{min})}, 2^{(l_{max})}]\f$. We do this by defining the parameters \f$\Sigma_{l} = \frac{1}{n^2} \sum_{i,j=1}^{n} K_{\epsilon_l}^{(ij)}\f$ and \f$\Sigma_{\epsilon_l}^{\prime} = (\log{(\Sigma_{\epsilon_{l+1}})}-\log{(\Sigma_{\epsilon_l})})/\log{(2)}\f$

  The optimal \f$\epsilon\f$ maximizes \f$\Sigma_{\epsilon}^{\prime}\f$ and the corresponding manifold dimension estimate is \f$m = 2 \max{(\Sigma_{\epsilon}^{\prime})}\f$.

  \return The density estimation at each sample \f$\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})\f$

  @param[in] tune <tt>true</tt> (default): Tune the bandwidth parameter values; <tt>false</tt>: use the stored parameters
  \return The density estimate \f$\psi^{(i)}\f$
  */
  Eigen::VectorXd Estimate(bool const tune = true);

  /// Estimate the density at each sample (given bandwidth parameter \f$r_i\f$)
  /**
  This function should be used if we have already computed the bandwidth parameter \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$.

  This function computes the kernel marix \f$\boldsymbol{K}_{\epsilon}\f$. We then estimate the density at each sample as
  \f{equation*}{
  \psi^{(i)}_{\epsilon} = \sum_{j=1}^{n} \frac{K_{\epsilon}^{(ij)}}{n (\pi \epsilon r_i^2)^{m/2}},
  \f}
  where \f$m\f$ is the manifold dimension.

  If <tt>tune</tt> is <tt>true</tt>, we must choose the optimal \f$\epsilon\f$ for the density estimation. Given the bandwidth exponent \f$l \in [l_{min}, l_{max}]\f$, we choose \f$\epsilon \in [2^{(l_{min})}, 2^{(l_{max})}]\f$. We do this by defining the parameters \f$\Sigma_{l} = \frac{1}{n^2} \sum_{i,j=1}^{n} K_{\epsilon_l}^{(ij)}\f$ and \f$\Sigma_{\epsilon_l}^{\prime} = (\log{(\Sigma_{\epsilon_{l+1}})}-\log{(\Sigma_{\epsilon_l})})/\log{(2)}\f$

  The optimal \f$\epsilon\f$ maximizes \f$\Sigma_{\epsilon}^{\prime}\f$ and the corresponding manifold dimension estimate is \f$m = 2 \max{(\Sigma_{\epsilon}^{\prime})}\f$.

  @param[in] squaredBandwidth The bandwidth squared parameter \f$r_i = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  \return The density estimation at each sample \f$\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})\f$

  @param[in] squaredBandwidth The squared bandwidth parameter \f$r_i^2\f$
  @param[in] tune <tt>true</tt> (default): Tune the bandwidth parameter values; <tt>false</tt>: use the stored Parameters
  \return The density estimate \f$\psi^{(i)}\f$
  */
  Eigen::VectorXd Estimate(Eigen::Ref<const Eigen::VectorXd> const& squaredBandwidth, bool const tune = true);

private:

  /// Tune the manifold dimension \f$m\f$?
  const bool tuneManifoldDimension;

  /// The default values for the spi::NumericalSolvers::DensityEstimation class.
  struct DefaultParameters {
    /// By default, do not tune the manifold dimension
    inline static const bool tuneManifoldDimension = false;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
