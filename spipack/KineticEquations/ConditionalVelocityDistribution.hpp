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
    @param[in] massDensity The mass density \f$\mu\f$
    @param[in] velocity The expected velocity \f$\boldsymbol{u}\f$
    @param[in] velocityDivergence The velocity divergence \f$\nabla_{\boldsymbol{x}} \cdot \boldsymbol{u}\f$
    @param[in] logMassDensityGrad The gradient of the log mass density \f$\nabla_{\boldsymbol{x}} \log{(\mu)}\f$
    */
    MacroscaleInformation(double const massDensity, Eigen::VectorXd const& velocity, double const velocityDivergence, Eigen::VectorXd const& logMassDensityGrad);

    virtual ~MacroscaleInformation() = default;

    /// The mass density \f$\mu\f$
    double massDensity;

    /// The expected velocity \f$\boldsymbol{u}\f$
    Eigen::VectorXd velocity;

    /// The velocity divergence \f$\nabla_{\boldsymbol{x}} \cdot \boldsymbol{u}\f$
    double velocityDivergence;

    // The gradient of the log mass density \f$\nabla_{\boldsymbol{x}} \log{(\mu)}\f$
    Eigen::VectorXd logMassDensityGrad;
  private:
  };

  /// Construct the conditional velocity distribution by sampling a random variable from \f$\psi\f$
  /**
  @param[in] macroLoc The location in the macroscale domain
  @param[in] rv The random variable that we wish to sample
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(Eigen::VectorXd const& macroLoc, std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  /// Construct the conditional velocity distribution given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] macroLoc The location in the macroscale domain
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(Eigen::VectorXd const& macroLoc, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  /// Construct the conditional velocity distribution given the samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] macroLoc The location in the macroscale domain
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(Eigen::VectorXd const& macroLoc, std::shared_ptr<spi::Tools::NearestNeighbors> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  virtual ~ConditionalVelocityDistribution() = default;

  /// Get the number of samples in the collection
  /**
  \return The number of samples.
  */
  std::size_t NumSamples() const;

  /// Get the state dimension of the samples
  /**
  \return The state dimension.
  */
  std::size_t StateDim() const;

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
  @param[in] saveInitialConditions <tt>true</tt> (default): save the initial conditions to file (only if we are writing to file); <tt>false</tt>: do not save the initial conditions (even if we are saving the other states to file). Useful if we are repeatedly calling this function using the previous end state as the intial state.
  */
  void Run(double const nextTime, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo, bool const saveInitialConditions = true);

  /// Get the \f$i^{th}\f$ point \f$\boldsymbol{v}^{(i)}\f$ from the point cloud.
  /**
  @param[in] i We want to get this point
  \return The \f$i^{th}\f$ point from the point cloud.
  */
  Eigen::Ref<Eigen::VectorXd const> Point(std::size_t const i) const;

  /// Get the \f$i^{th}\f$ point \f$\boldsymbol{v}^{(i)}\f$ from the point cloud.
  /**
  @param[in] i We want to get this point
  \return The \f$i^{th}\f$ point from the point cloud.
  */
  Eigen::Ref<Eigen::VectorXd> Point(std::size_t const i);

  /// Get the normalizing constant
  /**
  \return The normalizing constant
  */
  double NormalizingConstant() const;

  /// Compute the external acceleration
  /**
  By default, the external acceleration is zero.
  @param[in] vel The macro-scale velocity
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  virtual Eigen::VectorXd ExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& vel, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const;

  /// Compute the expected external acceleration
  /**
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  Eigen::VectorXd ExpectedExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const;

  /// The covariance given the expected velocity
  /**
  @param[in] expectedVel The expected velocity \f$\boldsymbol{u}\f$
  \return The covariance \f$\Sigma\f$
  */
  Eigen::MatrixXd Covariance(Eigen::Ref<const Eigen::VectorXd> const& mean) const;

  /// The expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$
  /**
  @param[in] expectedVel The expected velocity \f$\boldsymbol{u}\f$
  \return The expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$
  */
  double ExpectedEnergy(Eigen::Ref<const Eigen::VectorXd> const& expectedVel) const;

protected:

  /// The collision rate function \f$\ell\f$
  /**
  Compute the probability that a particle $\boldsymbol{v}^{(i)}\f$ collides. By default assume the collision rate function is constant \f$beta_0 = 1\f$.
  @param[in] x The macro-scale position \f$\boldsymbol{x}\f$
  @param[in] v The macro-scale velocity \f$\boldsymbol{v}\f$
  @param[in] t The macro-scale time \f$t\f$
  */
  virtual double CollisionRateFunction(Eigen::Ref<const Eigen::VectorXd> const& x, Eigen::Ref<const Eigen::VectorXd> const& v, double const t);

  /// Sample a vector from the unit hypersphere
  /**
  By default return \f$\boldsymbol{w} = (\boldsymbol{v}-\boldsymbol{v^{\prime}})/\|\boldsymbol{v}-\boldsymbol{v^{\prime}}\|\f$
  @param[in] v One of the velocity vectors \f$\boldsymbol{v}\f$
  @param[in] vprime The second velocity vector $\boldsymbol{v^{\prime}}\f$
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  @param[in] t The macro-scale time \f$t\f$
  \return A vector on the unit hypersphere \f$\boldsymbol{w}\f$
  */
  virtual Eigen::VectorXd SampleUnitHypersphere(Eigen::Ref<const Eigen::VectorXd> const& v, Eigen::Ref<const Eigen::VectorXd> const& vprime, Eigen::Ref<const Eigen::VectorXd> const& x, double const t) const;

  /// The post collision velocity function
  /**
  By default, assume perfectly elastic collisions.
  @param[in] v One of the velocity vectors \f$\boldsymbol{v}\f$
  @param[in] vprime The second velocity vector $\boldsymbol{v^{\prime}}\f$
  @param[in] w The angle vector between the two samples
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  @param[in] t The macro-scale time \f$t\f$
  \return A vector on the unit hypersphere \f$\boldsymbol{w}\f$
  */
  virtual double PostCollisionFunction(Eigen::Ref<const Eigen::VectorXd> const& v, Eigen::Ref<const Eigen::VectorXd> const& vprime, Eigen::Ref<const Eigen::VectorXd> const& w, Eigen::Ref<const Eigen::VectorXd> const& x, double const t) const;

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
    @param[in] dim The dimension of the state space
    */
    RescaledMacroscaleInformation(std::shared_ptr<const MacroscaleInformation> const& macroInfo, double const alpha, std::size_t const dim);

    virtual ~RescaledMacroscaleInformation() = default;

    /// The acceleration at each sample
    Eigen::MatrixXd acceleration;
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
  std::shared_ptr<ConditionalVelocityDistribution::RescaledMacroscaleInformation> InterpolateMacroscaleInformation(double const microT, std::shared_ptr<const ConditionalVelocityDistribution::MacroscaleInformation> const& nextMacroInfo);

  /// Update the normalizing constant
  /**
  @param[in] macroDelta The macro-scale timestep size
  @param[in] microDelta The micro-scale timestep size
  @param[in] prevInfo The macro-scale information at the previous micro-scale timestep
  @param[in] currInfo The macro-scale information at the current macro-scale timestep
  */
  void UpdateNormalizingConstant(double const macroDelta, double const microDelta, std::shared_ptr<const RescaledMacroscaleInformation> const& prevInfo, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo);

  /// Compute the collision step
  /**
  @param[in] macroDelta The macro-scale timestep size
  @param[in] microDelta The micro-scale timestep size
  @param[in] macroTime The macro-scale time
  @param[in] currInfo The macro-scale information at the current macro-scale timestep
  */
  void CollisionStep(double const macroDelta, double const microDelta, double const macroTime, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo);

  /// Update the particel velocities
  /**
  @param[in] macroDelta The macro-scale timestep size
  @param[in] microDelta The micro-scale timestep size
  @param[in] macroTime The macro-scale time
  @param[in] currInfo The macro-scale information at the current macro-scale timestep
  */
  void ConvectionStep(double const macroDelta, double const microDelta, double const macroTime, std::shared_ptr<RescaledMacroscaleInformation> const& currInfo);

  /// Update the particel velocities
  /**
  @param[in] macroDelta The macro-scale timestep size
  @param[in] macroTime The macroscale time
  @param[in] currInfo The macroscale information at the current macro-scale timestep
  */
  void ComputeAcceleration(double const macroDelta, double const macroTime, std::shared_ptr<RescaledMacroscaleInformation> const& currInfo);

  /// Solve the weighted Laplace problem to address source/sink term
  /**
   The effective acceleration due to the source/sink term at each point is stored in the rescaled macro-scale information.
  @param[in,out] currInfo The macroscale information at the current macro-scale timestep
  */
  void WeightedPoissonGradient(std::shared_ptr<RescaledMacroscaleInformation> const& currInfo) const;

  /// Write the samples to file
  /**
  @param[in] macroTime The macro-scale time
  @param[in] currInfo The macro-scale information at the current macro-scale timestep
  @param[in] file The file where we want to store the data.
  @param[in] dataset The data set where we put the data (defaults to <tt>"/"</tt>)
  */
  void WriteToFile(double const macroTime, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo, std::string const& dataset = "/") const;

  /// Compute the expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$ recursively
  /**
  @param[in] first The first sample index
  @param[in] last The last sample index
  @param[in] expectedVel The expected velocity \f$\boldsymbol{u}\f$
  \return The expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$
  */
  double ExpectedEnergy(std::size_t const first, std::size_t const last, Eigen::Ref<const Eigen::VectorXd> const& expectedVel) const;

  /// Compute the expected external acceleration
  /**
  @param[in] first The first sample index
  @param[in] last The last sample index
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  Eigen::VectorXd ExpectedExternalAcceleration(std::size_t const first, std::size_t const last, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const;

  /// The macroscale location \f$\boldsymbol{x}\f$.
  const Eigen::VectorXd macroLoc;

  /// The samples that we use to represent the distribution
  std::shared_ptr<spi::Tools::NearestNeighbors> samples;

  /// The Kolmogorov operator that allows us to compute augmented acceleration
  std::shared_ptr<spi::NumericalSolvers::KolmogorovOperator> kolmogorov;

  /// The macro-scale timestep
  double currentTime;

  /// The non-dimensional parameter \f$\varepsilon\f$
  const double varepsilon;

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

  /// The noise scaling for the external acceleration
  /**
  Add \f$a_0 \delta \boldsymbol{\tilde{A}}\f$, where \f$\boldsymbol{tilde{A}} \sim \mathcal{N}(\cdot; \boldsymbol{0}, \boldsymbol{I})\f$.
  */
  const double accelerationNoiseScale;

  /// The number of timesteps before we re-tune the Kolmogorov operator numerical parameters
  const std::size_t retuneFreq;

  /// The number of timesteps in the micro-scale simulation
  const std::size_t numTimesteps;

  /// The macro-scale information from the previous macro-scale timestep
  /**
  The information stored here is rescaled into the micro-scale coordinates.
  */
  std::shared_ptr<RescaledMacroscaleInformation> prevMacroInfo;

  /// The normalizing constant
  double normalizingConstant;

  /// The output filename
  const std::string filename;

  /// The default parameter values for spi::KineticEquations::ConditionalVelocityDistribution
  struct DefaultParameters {
    /// The default current time is \f$0\f$.
    inline static const double currentTime = 0.0;

    /// The default value for the non-dimensional parameter is \f$\varepsilon = \infty\f$.
    /**
    Setting the \f$\varepsilon\f$ parameter to \f$\infty\f$ "turns off" the collision step
    */
    inline static const double varepsilon = std::numeric_limits<double>::infinity();

    /// The default external acceleration rescaling is \f$1\f$.
    inline static const double alpha = 1.0;

    /// The default numerical timestep parameter is \f$\theta=0.5\f$
    inline static const double theta = 0.5;

    /// The default acceleration noise scale parameter is \f$1\f$.
    inline static const double accelerationNoiseScale = 1.0;

    /// The default number of timesteps before we retune the Kolmogorov operator parameters is a \f$10\f$.
    inline static const std::size_t retuneFreq = 10;

    /// The default number of timesteps is \f$10\f$
    inline static const std::size_t numTimesteps = 10;

    /// The default value for the initial normalizing constant is \f$1.0\f$.
    inline static const double initialNormalizingConstant = 1.0;

    /// The default output filename is empty
    inline static const std::string filename = "";
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace KineticEquations
} // namespace spi

#endif
