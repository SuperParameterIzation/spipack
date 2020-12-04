#ifndef CONDITIONALVELOCITYDISTRIBUTION_HPP_
#define CONDITIONALVELOCITYDISTRIBUTION_HPP_

#include "spipack/NumericalSolvers/SampleRepresentation/KolmogorovOperator.hpp"

namespace spi {
namespace KineticEquations{

class ConditionalVelocityDistribution {
public:

  /// The macro-scale information at this point (in the macro-scale coordinate system)
  class MacroscaleInformation {
  public:
    /// Construct the defaults
    MacroscaleInformation();

    /// Construct from existing data
    /**
    @param[in] velocityDivergence The velocity divergence \f$\nabla_{\boldsymbol{x}} \cdot \boldsymbol{u}\f$
    */
    MacroscaleInformation(double const velocityDivergence);

    virtual ~MacroscaleInformation() = default;

    /// The velocity divergence \f$\nabla_{\boldsymbol{x}} \cdot \boldsymbol{u}\f$
    /**
    Defaults to \f$0\f$.
    */
    double velocityDivergence = 0.0;
  private:
  };

  /// Construct the conditional velocity distribution by sampling a random variable from \f$\psi\f$
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  /// Construct the conditional velocity distribution given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  /// Construct the conditional velocity distribution given the samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(std::shared_ptr<spi::Tools::NearestNeighbors> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  virtual ~ConditionalVelocityDistribution() = default;

  /// Get the number of samples in the collection
  /**
  \return The number of samples.
  */
  std::size_t NumSamples() const;

  /// Get the current time
  /**
  \return The current time
  */
  double CurrentTime() const;

  /// Get the numerical timestep parameter
  /**
  \return The numerical timestep parameter
  */
  double TimestepParameter() const;

  /// Get the number of timesteps
  /**
  \return The number of timesteps
  */
  std::size_t NumTimesteps() const;

  /// Get the external acceleration rescaling (\f$\alpha\f$ parameter)
  /**
  \return The external acceleration rescaling \f$\alpha\f$
  */
  double ExternalAccelerationRescaling() const;

  /// Get the current macro-scale information
  /**
  This function rescales the stored information into macro-scale coordinates.
  \return The macro-scale information in macro-scale coordinates.
  */
  std::shared_ptr<const ConditionalVelocityDistribution::MacroscaleInformation> MacroscaleInfo() const;

  /// Run the model for a macro-scale timestep
  /**
  @param[in] nextTime We need to update the samples from the current time to this time
  @param[in] nextMacroInfo The macro-scale information at the updated macro-scale timestep (in the macro-scale coordinates)
  */
  void Run(double const nextTime, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo);

  /// Get the \f$i^{th}\f$ point \f$\boldsymbol{v}^{(i)}\f$ from the point cloud.
  /**
  @param[in] i We want to get this point
  \return The \f$i^{th}\f$ point from the point cloud.
  */
  Eigen::Ref<Eigen::VectorXd const> Point(std::size_t const i) const;

  /// Get the normalizing constant
  /**
  \return The normalizing constant
  */
  double NormalizingConstant() const;

private:

  /// The macro-scale information at this point (in the micro-scale coordinate system)
  class RescaledMacroscaleInformation : public MacroscaleInformation {
  public:
    /// Construct from defaults
    RescaledMacroscaleInformation();

    /// Construct from existing data
    /**
    @param[in] macroInfo The macro-scale information in the macro-scale coordinate system
    @param[in] alpha The external acceleration rescaling parameter
    */
    RescaledMacroscaleInformation(std::shared_ptr<const MacroscaleInformation> const& macroInfo, double const alpha);

    virtual ~RescaledMacroscaleInformation() = default;
  private:
  };

  /// Update the options for the Kolmogorov options to make sure all the options are valid
  /**
  @param[in] options The user-input options
  \return The updated options that ensure validity
  */
  static YAML::Node KolmogorovOptions(YAML::Node options, std::size_t const dim);

  /// Compute the current macroscale information
  /**
  Given the micro-scale time, interpolate between the previous macro-scale information and the next macro-scale information
  @param[in] microT The micro-scale time
  @param[in] nextMacroInfo The macro-scale information at the next macro-scale timestep
  \return The macro-scale information interpolated onto the micro-scale time
  */
  std::shared_ptr<const ConditionalVelocityDistribution::RescaledMacroscaleInformation> InterpolateMacroscaleInformation(double const microT, std::shared_ptr<const ConditionalVelocityDistribution::MacroscaleInformation> const& nextMacroInfo);

  /// Update the normalizing constant
  /**
  @param[in] macroDelta The macro-scale timestep size
  @param[in] microDelta The micro-scale timestep size
  @param[in] prevInfo The macro-scale information at the previous micro-scale timestep
  @param[in] currInfo The macro-scale information at the current macro-scale timestep
  */
  void UpdateNormalizingConstant(double const macroDelta, double const microDelta, std::shared_ptr<const RescaledMacroscaleInformation> const& prevInfo, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo);

  /// The samples that we use to represent the distribution
  std::shared_ptr<spi::Tools::NearestNeighbors> samples;

  /// The Kolmogorov operator that allows us to compute augmented acceleration
  std::shared_ptr<spi::NumericalSolvers::KolmogorovOperator> kolmogorov;

  /// The macro-scale timestep
  double currentTime;

  /// The rescaling for the acceleration
  /**
  The rescaling \f$\alpha\f$ also defines the spatial and velocity rescaling \f$\zeta = \alpha\f$ and \f$\xi = \alpha\f$.
  */
  const double alpha;

  /// The numerical timestep parameter
  /**
  \f$\theta=0\f$ corresponds to a fully implicit timestepping algorithm; we require that \f$\theta \in [0,1]\f$.
  */
  const double theta;

  /// The number of timesteps in the micro-scale simulation
  const std::size_t numTimesteps;

  /// The macro-scale information from the previous macro-scale timestep
  /**
  The information stored here is rescaled into the micro-scale coordinates.
  */
  std::shared_ptr<const RescaledMacroscaleInformation> prevMacroInfo;

  /// The normalizing constant
  double normalizingConstant;

  /// The default parameter values for spi::KineticEquations::ConditionalVelocityDistribution
  struct DefaultParameters {
    /// The default current time is \f$0\f$.
    inline static const double currentTime = 0.0;

    /// The default external acceleration rescaling is \f$1\f$.
    inline static const double alpha = 1.0;

    /// The default numerical timestep parameter is \f$\theta=0.5\f$
    inline static const double theta = 0.5;

    /// The default number of timesteps is \f$10\f$
    inline static const std::size_t numTimesteps = 10;

    /// The default value for the initial normalizing constant is \f$1.0\f$.
    inline static const double initialNormalizingConstant = 1.0;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace KineticEquations
} // namespace spi

#endif
