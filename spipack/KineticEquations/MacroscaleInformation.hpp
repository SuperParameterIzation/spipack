#ifndef MACROSCALEINFORMATION_HPP_
#define MACROSCALEINFORMATION_HPP_

#include <cstddef>

#include <Eigen/Core>

namespace spi {
namespace KineticEquations {

/// Is the information stored in macro- or micro-scale coordinates?
enum Coordinates {
  /// The samples are stored in macro-scale coordinates.
  MACRO,

  /// The samples are stored in micro-scale coordinates.
  MICRO
};

class MacroscaleInformation {
public:

  /// The spatial dimension (is always \f$2\f$).
  inline constexpr static std::size_t dim = 2;

  /// Construct the macro-scale information with default values
  MacroscaleInformation();

  /// Construct the macro-scale information with given values
  /**
  @param[in] massDensity The mass density
  @param[in] logMassDensityGrad The gradient of the log mass density
  @param[in] velocity The velocity
  @param[in] velocityDiv The divergence of the velocity
  @param[in] coordinates The coordinate system that the data is stored in
  */
  MacroscaleInformation(double const massDensity, Eigen::Matrix<double, dim, 1> const& logMassDensityGrad, Eigen::Matrix<double, dim, 1> const& velocity, double const velocityDiv, Coordinates const& coordinates);

  virtual ~MacroscaleInformation() = default;

  /// The mass density \f$\mu(\boldsymbol{x}; t)\f$ (macro) or \f$\phi(\boldsymbol{X}; T)\f$ (micro).
  /**
  Defaults to \f$\mu(\boldsymbol{x}; t) = 1\f$.
  */
  double massDensity = 1.0;

  /// The gradient of the log-density \f$\nabla_{\boldsymbol{x}} \log{(\mu(\boldsymbol{x}; t))}\f$ (macro) or \f$\nabla_{\boldsymbol{X}} \log{(\phi(\boldsymbol{X}; T))}\f$ (micro)
  /**
  Defaults to \f$\boldsymbol{0}\f$.
  */
  Eigen::Matrix<double, dim, 1> logMassDensityGrad = Eigen::Matrix<double, dim, 1>::Zero(dim, 1);

  /// The expected velocity \f$\boldsymbol{u}\f$ (macro) or \f$\boldsymbol{U}\f$ (micro)
  /**
  Defaults to \f$\boldsymbol{0}\f$.
  */
  Eigen::Matrix<double, dim, 1> velocity = Eigen::Matrix<double, dim, 1>::Zero(dim, 1);

  /// The divergence of the velocity
  /**
  Defaults to \f$0\f$.
  */
  double velocityDiv = 0.0;

  /// Which coordinate system is this information stored in?
  /**
  Defaults to the macro-scale.
  */
  Coordinates coordinates = Coordinates::MACRO;

  // This method lets cereal know which data members to serialize
  template<class Archive>
  void serialize(Archive & archive) { archive(massDensity, logMassDensityGrad, velocity, velocityDiv, coordinates); }

private:

};

} // namespace KineticEquations
} // namespace spi

#endif
