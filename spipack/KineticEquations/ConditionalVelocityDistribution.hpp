#ifndef CONDITIONALVELOCITYDISTRIBUTION_HPP_
#define CONDITIONALVELOCITYDISTRIBUTION_HPP_

#include <Sacado.hpp>

#include <MUQ/SamplingAlgorithms/SampleGraphs/KolmogorovOperator.h>

#include <boost/property_tree/ptree.hpp>

#include <Eigen/Sparse>

#include "spipack/Tools/NearestNeighbors.hpp"

#include "spipack/KineticEquations/MacroscaleInformation.hpp"

namespace spi {
namespace KineticEquations{

/// Forward declaration of the right hand side for the convection step
class ConvectionRHS;

/// The micro-scale model for the rescaled conditional velocity distribution \f$\psi(\boldsymbol{V} \vert \boldsymbol{X}; T) = \rho \hat{\psi}(\boldsymbol{V} \vert \boldsymbol{X}; T)\f$
/**
\f{eqnarray*}{
\partial_T \rho &=& \delta \rho \nabla_{\boldsymbol{X}} \cdot \boldsymbol{U} \\
\partial_T \hat{\psi} + \delta \nabla_{\boldsymbol{V}} \cdot (\hat{\psi} \boldsymbol{A}) &=& \phi (\rho^{-1} Q[\rho \hat{\psi}] - L\hat{\psi}) - \delta (\boldsymbol{V} + \boldsymbol{\hat{U}} - \boldsymbol{U}) \cdot \nabla _{\boldsymbol{X}} \log{(\phi)}
\f}
<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"NearestNeighbors"   | <tt>YAML::Node</tt> | - | Options for the spi::Tools::NearestNeighbors object.   |
"KolmogorovOptions"   | <tt>YAML::Node</tt> | - | Options for the the spi::NumericalSolvers::KolmogorovOperator object. |
"CurrentTime" | <tt>double</tt> | <tt>0.0</tt> | The macro-scale time when this object is constructed
"NondimensionalParameter" | <tt>double</tt> | \f$\infty\f$ | The nondimensional parameter that defines the collision rate scaling
"AccelerationNoiseScale" | <tt>double</tt> | <tt>1.0</tt> | The amount of noise we add to the external acceleration---scaled also by the macro-scale timestep \f$\delta\f$.
"KolmogorovRetuneFrequency" | <tt>std::size_t</tt> | \f$\infty\f$ | The number of timesteps before we retune the numerical parameters for the spi::NumericalSolvers::KolmogorovOperator.
"NumTimesteps" | <tt>std::size_t</tt> | <tt>10</tt> | The number of timesteps for the micro-scale model
"OutputFilename" | <tt>std::string</tt> | <tt>""</tt> | The file where we output the data. If this string is empty, do not output any data.
*/
class ConditionalVelocityDistribution : public std::enable_shared_from_this<ConditionalVelocityDistribution> {
public:
  /// The right hand side for the convection step is a friend of this class
  friend ConvectionRHS;

  /// Construct the conditional velocity distribution by sampling a random variable from \f$\psi\f$
  /**
  @param[in] macroLoc The location in the macroscale domain
  @param[in] rv The random variable that we wish to sample
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& macroLoc, std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  /// Construct the conditional velocity distribution given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] macroLoc The location in the macroscale domain
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& macroLoc, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

  /// Construct the conditional velocity distribution given the samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] macroLoc The location in the macroscale domain
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  ConditionalVelocityDistribution(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& macroLoc, std::shared_ptr<spi::Tools::NearestNeighbors> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options);

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

  /// Get the current macro-scale information
  /**
  This function rescales the stored information into macro-scale coordinates.
  \return The macro-scale information in macro-scale coordinates.
  */
  std::shared_ptr<const MacroscaleInformation> MacroscaleInfo() const;

  /// Run the model for a macro-scale timestep
  /**
  @param[in] nextTime We need to update the samples from the current time to this time
  @param[in] nextMacroInfo The macro-scale information at the updated macro-scale timestep (in the macro-scale coordinates)
  @param[in] saveInitialConditions <tt>true</tt> (default): save the initial conditions to file (only if we are writing to file); <tt>false</tt>: do not save the initial conditions (even if we are saving the other states to file). Useful if we are repeatedly calling this function using the previous end state as the intial state.
  */
  void Run(double const nextTime, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo);

  /// Get the \f$i^{th}\f$ point \f$\boldsymbol{v}^{(i)}\f$ from the point cloud.
  /**
  @param[in] i We want to get this point
  \return The \f$i^{th}\f$ point from the point cloud.
  */
  Eigen::Ref<Eigen::Matrix<double, MacroscaleInformation::dim, 1> const> Point(std::size_t const i) const;

  /// Get the \f$i^{th}\f$ point \f$\boldsymbol{v}^{(i)}\f$ from the point cloud.
  /**
  @param[in] i We want to get this point
  \return The \f$i^{th}\f$ point from the point cloud.
  */
  Eigen::Ref<Eigen::Matrix<double, MacroscaleInformation::dim, 1> > Point(std::size_t const i);

  /// Compute the external acceleration
  /**
  By default, the external acceleration is zero.
  @param[in] vel The macro-scale velocity
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  virtual Eigen::Matrix<double, MacroscaleInformation::dim, 1> ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vel, double const time) const;

  /// Compute the expected external acceleration
  /**
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> ExpectedExternalAcceleration(double const time);

  /// Get the covariance matrix of the current samples
  /**
  Computes and stores the covarinace, if it has already been computed then just return the stored covariance
  \return The covariance matrix
  */
  Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> Covariance();

  /// Get the skew vector of the current samples
  /**
  Computes and stores the skew, if it has already been computed then just return the stored skew
  \return The skew vector
  */
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> Skew();

  /// Compute the skew vector of the current samples
  /**
  Recursively computes the skew
  @param[in] first The index of the first sample used to compute the skew
  @param[in] last The index of the last sample used to compute the skew
  \return The skew vector
  */
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> Skew(std::size_t const first, std::size_t const last) const;

  /// The expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$
  /**
  @param[in] expectedVel The expected velocity \f$\boldsymbol{u}\f$
  \return The expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$
  */
  double ExpectedEnergy(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& expectedVel);

  /// The macro-scale location
  /**
  \return The macro-scale location
  */
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> MacroscaleLocation() const;

  /// Compute the sample mean
  /**
  This should the macro-scale velocity, so this function is really only for diagonstics and testing
  \return The sample mean \f$n^{-1} \sum_{i=1}^{n} \boldsymbol{v}^{(i)}\f$.
  */
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> SampleMean() const;

protected:

  /// The collision rate function \f$\ell\f$
  /**
  Compute the probability that a particle $\boldsymbol{v}^{(i)}\f$ collides. By default assume the collision rate function is constant \f$beta_0 = 1\f$.
  @param[in] vi The macro-scale velocity \f$\boldsymbol{v}^{(i)}\f$
  @param[in] t The macro-scale time \f$t\f$
  */
  virtual double CollisionRateFunction(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& vi, double const t) const;

  /// Sample a vector from the unit hypersphere
  /**
  By default return \f$\boldsymbol{w} = (\boldsymbol{v}-\boldsymbol{v^{\prime}})/\|\boldsymbol{v}-\boldsymbol{v^{\prime}}\|\f$
  @param[in] v One of the velocity vectors \f$\boldsymbol{v}\f$
  @param[in] vprime The second velocity vector $\boldsymbol{v^{\prime}}\f$
  @param[in] t The macro-scale time \f$t\f$
  \return A vector on the unit hypersphere \f$\boldsymbol{w}\f$
  */
  virtual Eigen::Matrix<double, MacroscaleInformation::dim, 1> SampleUnitHypersphere(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& v, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vprime, double const t) const;

  /// The post collision velocity function
  /**
  By default, assume perfectly elastic collisions.
  @param[in] v One of the velocity vectors \f$\boldsymbol{v}\f$
  @param[in] vprime The second velocity vector $\boldsymbol{v^{\prime}}\f$
  @param[in] w The angle vector between the two samples
  @param[in] t The macro-scale time \f$t\f$
  \return A vector on the unit hypersphere \f$\boldsymbol{w}\f$
  */
  virtual double PostCollisionFunction(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& v, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vprime, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& w, double const t) const;

  /// The macroscale location \f$\boldsymbol{x}\f$.
  const Eigen::Matrix<double, MacroscaleInformation::dim, 1> macroLoc;

private:

  /// Create the Kolmogorov operator
  /**
  @param[in] options The options to construct the entire ConditionalVelocityDistribution
  */
  void InitializeKolmogorovOperator(YAML::Node const& options);

  /// Shift the samples so they match the expected velocity from the macro-scale information
  /**
  @param[in] macroInfo Shift the samples so the mean matches this velocity (this is in rescaled coordinates)
  */
  void ShiftSamples(std::shared_ptr<const MacroscaleInformation> const& macroInfo);

  void CollisionStep(double const delta, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo);

  /// Compute the expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$ recursively
  /**
  @param[in] first The first sample index
  @param[in] last The last sample index
  @param[in] expectedVel The expected velocity \f$\boldsymbol{u}\f$
  \return The expected energy \f$\int_{\mathcal{V}} \frac{1}{2} \boldsymbol{v} \cdot \boldsymbol{v} \psi \, d\boldsymbol{v}\f$
  */
  double ExpectedEnergy(std::size_t const first, std::size_t const last, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& expectedVel) const;

  /// Compute the expected external acceleration
  /**
  @param[in] first The first sample index
  @param[in] last The last sample index
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> ExpectedExternalAcceleration(std::size_t const first, std::size_t const last, double const time) const;

  /// Update the samples (both from the distribution and the tilted plane)
  /**
  @param[in] currInfo The macroscale information at the current macro-scale timestep
  @param[in] density An estimate of the density at each sample
  @param[in] stepsize The timestep size
  @param[in,out] tiltedSamples The samples that move off the \f$\boldsymbol{x} = \boldsymbol{\hat{x}}\f$ plane
  */
  void UpdateSamples(std::shared_ptr<const MacroscaleInformation> const& currInfo, Eigen::VectorXd const& density, double const stepsize, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& tiltedSamples);

  /// Evaluated the tilted conditional density at the point parameterized by the samples
  /**
  @param[in] currInfo The macroscale information at the current macro-scale timestep
  @param[in] density An estimate of the density at each sample
  @param[in] stepsize The timestep size
  @param[in] tiltedSamples The samples that move off the \f$\boldsymbol{x} = \boldsymbol{\hat{x}}\f$ plane
  \return An estimate of the tilted density using Nystrom's method
  */
  Eigen::VectorXd NystromMethod(std::shared_ptr<const MacroscaleInformation> const& currInfo, Eigen::VectorXd const& density, double const stepsize, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& tiltedSamples) const;

  /// Initialize the directional derivative \f$ \boldsymbol{v} \cdot \nabla_{\boldsymbol{x}} \nu \f$.
  void InitializeDirectionalDerivative();

  /// Evaluate the directional derivative at a point in the velocity space
  /**
  Defaults to returning 0.
  @param[in] velocity The velocity
  \return The directional derivative \f$ \boldsymbol{v} \cdot \nabla_{\boldsymbol{x}} \nu \f$.
  */
  virtual double DirectionalDerivative(Eigen::VectorXd const& velocity) const;

  /// The expected acceleration \f$\int \boldsymbol{A} \psi \, d \boldsymbol{V}\f$
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> expectedAcceleration;

  /// <tt>true</tt>: we have updated the expected acceleration to the current time; <tt>false</tt>: we need to recompute the expected acceleration
  bool computedExpectedAcceleration = false;

  /// The covariance of the samples
  Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> covariance;

  /// <tt>true</tt>: we have updated the covariance to the current time; <tt>false</tt>: we need to recompute the covariance
  bool computedCovariance = false;

  /// The skew of the samples
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> skew;

  /// <tt>true</tt>: we have updated the skew to the current time; <tt>false</tt>: we need to recompute the skew
  bool computedSkew = false;

  /// The expected energy
  double expectedEnergy;

  /// tt>true</tt>: we have updated the expected energy to the current time; <tt>false</tt>: we need to recompute the expected energy
  bool computedEnergy = false;

  /// The samples that we use to represent the distribution
  std::shared_ptr<spi::Tools::NearestNeighbors> samples;

  /// The Kolmogorov operator that allows us to compute augmented acceleration
  std::shared_ptr<muq::SamplingAlgorithms::KolmogorovOperator> kolmogorov;

  /// The directional derivative \f$ \boldsymbol{v} \cdot \nabla_{\boldsymbol{x}} \nu \f$.
  Eigen::VectorXd directionalDerivative;

  /// The macro-scale timestep
  double currentTime;

  /// The non-dimensional parameter \f$\varepsilon\f$
  const double varepsilon;

  /// The macro-scale information from the previous macro-scale timestep
  /**
  The information stored here is rescaled into the micro-scale coordinates.
  */
  std::shared_ptr<const MacroscaleInformation> prevMacroInfo;

  /// Options for the density estimation
  boost::property_tree::ptree densityEstimationOptions;

  struct DefaultParameters {
    /// The default current time is \f$0\f$.
    inline static const double currentTime = 0.0;

    /// The default value for the non-dimensional parameter is \f$\varepsilon = \infty\f$.
    /**
    Setting the \f$\varepsilon\f$ parameter to \f$\infty\f$ "turns off" the collision step
    */
    inline static const double varepsilon = std::numeric_limits<double>::infinity();
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace KineticEquations
} // namespace spi

#endif
