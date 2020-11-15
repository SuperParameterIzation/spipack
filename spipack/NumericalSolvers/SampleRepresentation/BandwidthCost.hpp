#ifndef BANDWIDTHCOST_HPP_
#define BANDWIDTHCOST_HPP_

#include <MUQ/Optimization/CostFunction.h>

#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

namespace spi {
namespace NumericalSolvers {

/// The cost function used to tune the bandwidth parameter \f$\epsilon\f$.
/**
If <tt>tune</tt> is <tt>true</tt>, we must choose the optimal \f$\epsilon\f$ for the density estimation. Given the bandwidth exponent \f$l \in [l_{min}, l_{max}]\f$, we choose \f$\epsilon \in [2^{(l_{min})}, 2^{(l_{max})}]\f$. We do this by defining the parameters \f$\Sigma_{l} = \frac{1}{n^2} \sum_{i,j=1}^{n} K_{\epsilon_l}^{(ij)}\f$ and \f$\Sigma_{\epsilon_l}^{\prime} = (\log{(\Sigma_{\epsilon_{l+1}})}-\log{(\Sigma_{\epsilon_l})})/\log{(2)}\f$

The optimal \f$\epsilon\f$ maximizes \f$\Sigma_{\epsilon}^{\prime}\f$ and the corresponding manifold dimension estimate is \f$m = 2 \max{(\Sigma_{\epsilon}^{\prime})}\f$.
*/
class BandwidthCost : public muq::Optimization::CostFunction {
public:

  /// Construct the cost function
  /**
  @param[in] bandwidth The vector \f$\boldsymbol{r}\f$ that defines the kernel matrix
  @param[in] samples A discrete sample representation of the density \f$\psi\f$
  */
  BandwidthCost(Eigen::Ref<const Eigen::VectorXd> const& bandwidth, std::shared_ptr<const SampleRepresentation> const& samples);

  virtual ~BandwidthCost() = default;
private:
  /// An implementation of the cost function \f$\Sigma_l^{\prime}\f$
  /**
  @param[in] input Only one input vector: a vector of length \f$1\f$ containing the candidate parameter \f$\epsilon\f$
  \return The cost function \f$\Sigma_l^{\prime}\f$
  */
  virtual double CostImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& input) override;

  /// An implementation of the cost function derivative \f$\frac{d}{d \epsilon}\Sigma_l^{\prime}\f$
  /**
  @param[in] inputDimWrt Since there is only one input vector, this must be \f$0\f$.
  @param[in] input Only one input vector: a vector of length \f$1\f$ containing the candidate parameter \f$\epsilon\f$
  */
  virtual void GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) override;

  /// The bandwidth vector \f$\boldsymbol{r}\f$ that defines the kernel matrix
  const Eigen::VectorXd bandwidth;

  /// A discrete sample representation of the density \f$\psi\f$
  std::shared_ptr<const SampleRepresentation> samples;

  /// The value of the bandwidth parameter \f$\epsilon_l = 2^{l}\f$ at the most recent function call.
  /**
  We store the results (because they are expensive to compute).
  */
  double eps = std::numeric_limits<double>::quiet_NaN();

  /// A cached value of the sum of the entries in the kernel matrix at \f$\epsilon_l = 2^{l}\f$
  double ksuml = std::numeric_limits<double>::quiet_NaN();

 /// A cached value of the sum of the entries in the kernel matrix at \f$\epsilon_{l+1} = 2^{l+1}\f$
 double ksumlp1 = std::numeric_limits<double>::quiet_NaN();
};

} // namespace NumericalSolvers
} // namespace spi

#endif
