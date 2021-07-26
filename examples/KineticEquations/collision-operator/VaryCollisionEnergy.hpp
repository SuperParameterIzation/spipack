#ifndef VARYCOLLISIONENERGY_HPP_
#define VARYCOLLISIONENERGY_HPP_

#include <MUQ/Utilities/RandomGenerator.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <spipack/KineticEquations/ConditionalVelocityDistribution.hpp>

class VaryCollisionEnergy : public spi::KineticEquations::ConditionalVelocityDistribution {
public:
  /// The dimension of the state space
  inline static constexpr std::size_t dim = 2;

  /// Construct the conditional velocity distribution with a linear external forcing
  /**
  @param[in] gamma The fraction of remaining kinetic energy after collision
  @param[in] options Setup options
  */
  inline VaryCollisionEnergy(double const gamma, YAML::Node const& options) : ConditionalVelocityDistribution(
    Eigen::Matrix<double, dim, 1>::Zero(dim),
    std::make_shared<muq::Modeling::Gaussian>(dim)->AsVariable(), DefaultMacroscaleInformation(),
    options),
    gamma(gamma)
    {}

  virtual ~VaryCollisionEnergy() = default;

  /**
  \return The macro-scale information with constat mass density, zero velocity, and zero derivatives.
  */
  inline static std::shared_ptr<MacroscaleInformation> DefaultMacroscaleInformation() {
    return std::make_shared<MacroscaleInformation>(1.0, Eigen::Matrix<double, dim, 1>::Zero(dim, 1), 0.0, Eigen::Matrix<double, dim, 1>::Zero(dim, 1));
  }

  /**
  \return The options for this simulation
  */
  inline static YAML::Node Options(std::string const& filename) {
    // the number of samples
    static const std::size_t n = 1000;

    // the number of timesteps
    static const std::size_t numTimesteps = 100;

    // options for the nearest neighbor search
    static YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;
    nnOptions["NumThreads"] = omp_get_max_threads();

    // set the options for the conditional velocity distribution
    YAML::Node options;
    options["NearestNeighbors"] = nnOptions;
    options["NumTimesteps"] = numTimesteps;
    options["OutputFilename"] = filename;
    options["AccelerationNoiseScale"] = 0.0;
    options["NondimensionalParameter"] = 1.0e-1;
    return options;
  }
private:
  /// Sample a vector from the unit hypersphere
  /**
  @param[in] v One of the velocity vectors \f$\boldsymbol{v}\f$
  @param[in] vprime The second velocity vector $\boldsymbol{v^{\prime}}\f$
  @param[in] t The macro-scale time \f$t\f$
  \return A vector on the unit hypersphere \f$\boldsymbol{w}\f$
  */
  inline virtual Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> SampleUnitHypersphere(Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& v, Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& vprime, double const t) const override {
    const double theta = 2.0*M_PI*muq::Utilities::RandomGenerator::GetUniform();
    Eigen::Matrix<double, dim, 1> w(dim, 1);
    w(0) = std::cos(theta); w(1) = std::sin(theta);
    return w;
  }

  /// The post collision velocity function
  /**
  @param[in] v One of the velocity vectors \f$\boldsymbol{v}\f$
  @param[in] vprime The second velocity vector $\boldsymbol{v^{\prime}}\f$
  @param[in] w The angle vector between the two samples
  @param[in] t The macro-scale time \f$t\f$
  \return A vector on the unit hypersphere \f$\boldsymbol{w}\f$
  */
  inline virtual double PostCollisionFunction(Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& v, Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& vprime, Eigen::Ref<const Eigen::Matrix<double, ConditionalVelocityDistribution::dim, 1> > const& w, double const t) const override {
    // the energy
    const double etimes4 = 2.0*(v.dot(v)+vprime.dot(vprime));

    // the perfectly elastic case
    const double sigmae = w.dot(vprime-v);
    const double sigmae2 = sigmae*sigmae;

    // the post collision velocity function
    return 0.5*(sigmae + std::copysign(std::sqrt(std::max(0.0, sigmae2-(1.0-gamma)*etimes4)), sigmae));
  }

  const double gamma = 1.0;
};

#endif
