#ifndef KINETICMODELS_HPP_
#define KINETICMODELS_HPP_

#include <Sacado.hpp>

#include <MUQ/Approximation/Regression/Regression.h>

#include <MUQ/SamplingAlgorithms/SampleCollection.h>

#include "spipack/Tools/NearestNeighbors.hpp"

#include "spipack/KineticEquations/ConditionalVelocityDistribution.hpp"

namespace spi {
namespace KineticEquations {

class KineticModels {
public:
  typedef Sacado::Fad::DFad<double> ADScalar;

  /// Get the kinetic models within a macro-scale simulation
  /**
  Get's the kinetic models stored within the <tt>ExaHyPE</tt> simulation.
  */
  static std::shared_ptr<KineticModels> ExaHyPEKineticModels();

  KineticModels();

  KineticModels(YAML::Node const& options);

  virtual ~KineticModels() = default;

  /**
  @param[in] loc The macro-scale location
  @param[in] mass The macro-scale mass density
  @param[in] velocity The macro-scale (expected) velocity
  @param[in] energy The macro-scale expected energy
  \return The flux matrix for the conservation equations
  */
  template<typename scalar>
  inline Eigen::Matrix<scalar, MacroscaleInformation::dim, 4> FluxMatrixImplmentation(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, scalar const mass, Eigen::Ref<const Eigen::Matrix<scalar, MacroscaleInformation::dim, 1> > const& velocity, scalar const energy) const {
    /*Eigen::Matrix<double, MacroscaleInformation::dim, 1> massGradient;
    Eigen::Matrix<double, MacroscaleInformation::dim, 1> pressureGradient;
    Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> velJac;
    if( !fields ) {
      Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 1> locAD;
      for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { locAD(i) = ADScalar(MacroscaleInformation::dim, i, loc(i)); }

      ADScalar massDensityAD = this->InitialMassDensity<ADScalar>(locAD);

      Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 1> velAD;
      Eigen::Matrix<double, MacroscaleInformation::dim, 1> vel;
      this->InitialExpectedVelocity<ADScalar>(locAD, velAD);

      Eigen::Matrix<double, MacroscaleInformation::dim, 1> energyGradient;
      ADScalar energyAD = this->InitialExpectedEnergy<ADScalar>(locAD);

      for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
        massGradient(i) = massDensityAD.dx(i);
        energyGradient(i) = energyAD.dx(i);
        vel(i) = velAD(i).val();
        for( std::size_t j=0; j<MacroscaleInformation::dim; ++j ) {
          velJac(i, j) = velAD(i).dx(j);
        }
      }
      pressureGradient = energyGradient - velJac*vel;
    } else {
      constexpr std::size_t nNeighs = 25;
      std::vector<std::pair<std::size_t, double> > neighbors;
      fields->FindNeighbors(loc, nNeighs, neighbors);
      assert(neighbors.size()==nNeighs);

      // the Vandermonde matrix and data for linear regression
      Eigen::Matrix<double, nNeighs, MacroscaleInformation::dim+1> vand;
      Eigen::Matrix<double, nNeighs, MacroscaleInformation::dim+2> data;
      // loop through the neighbors
      for( std::size_t i=0; i<nNeighs; ++i ) {
        auto sample = fields->GetSamplingState(neighbors[i].first);

        vand(i, 0) = 1.0;
        vand.row(i).tail(MacroscaleInformation::dim) = sample->state[0];

        // the mass
        data(i, 0) = std::max(0.0, sample->state[1] (0));

        // the velocity
        data.row(i).segment(1, MacroscaleInformation::dim) = sample->state[2];

        // the energy
        data(i, MacroscaleInformation::dim+1) = std::max(0.0, sample->state[3] (0));
      }

      // solve the linear regression
      Eigen::ColPivHouseholderQR<Eigen::Matrix<double, MacroscaleInformation::dim+1, MacroscaleInformation::dim+1> > qrfac;
      qrfac.compute(vand.transpose()*vand);
      assert(qrfac.isInvertible());
      const Eigen::Matrix<double, MacroscaleInformation::dim+1, MacroscaleInformation::dim+2> coeff = qrfac.solve(vand.transpose()*data);

      massGradient = coeff.col(0).tail(MacroscaleInformation::dim);

      //std::cout << mass << " " << coeff(0, 0) + massGradient.dot(loc) << std::endl;

      Eigen::Matrix<double, MacroscaleInformation::dim, 1> energyGradient = coeff.col(MacroscaleInformation::dim+1).tail(MacroscaleInformation::dim);

      Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> velJac = coeff.block(1, 1, MacroscaleInformation::dim, MacroscaleInformation::dim);

      Eigen::Matrix<double, MacroscaleInformation::dim, 1> vel = coeff.row(0).segment(1, MacroscaleInformation::dim);
      vel = vel + velJac.transpose()*loc;

      pressureGradient = energyGradient - velJac*vel;
    }

    Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> pressureDirs;
    const double pressureGradNorm = pressureGradient.norm();
    if( pressureGradNorm>1.0e-6 ) {
      //pressureGradient /= pressureGradNorm;
      pressureDirs(0, 0) = pressureGradient(1);
      pressureDirs(1, 0) = -pressureGradient(0);
      pressureDirs(0, 1) = pressureGradient(0);
      pressureDirs(1, 1) = pressureGradient(1);
    } else {
      pressureDirs = Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>::Identity(MacroscaleInformation::dim, MacroscaleInformation::dim);
    }
    //pressureDirs = Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>::Identity(MacroscaleInformation::dim, MacroscaleInformation::dim);

    Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> massDirs;
    const double massGradNorm = massGradient.norm();
    if( massGradNorm>1.0e-6 ) {
      massGradient /= massGradNorm;
      massDirs(0, 0) = massGradient(1);
      massDirs(1, 0) = -massGradient(0);
      massDirs(0, 1) = massGradient(0);
      massDirs(1, 1) = massGradient(1);
    } else {
      massDirs = Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>::Identity(MacroscaleInformation::dim, MacroscaleInformation::dim);
    }*/

    const scalar pressure = std::max(0.0, energy - 0.5*velocity.dot(velocity));

    //const Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = pressure*Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>::Identity(MacroscaleInformation::dim, MacroscaleInformation::dim);
    /*Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>::Zero(MacroscaleInformation::dim, MacroscaleInformation::dim);
    const double gamma = 0.9;
    cov(0, 0) = 2.0*gamma*pressure;
    cov(1, 1) = 2.0*(1.0-gamma)*pressure;*/

    /*std::cout << "mag: " << massGradientNorm << std::endl;
    std::cout << "proj mag: " << massGradient(0)*pressureDirs(0, 0) + massGradient(1)*pressureDirs(1, 0) << std::endl;
    std::cout << "gamma: " << (massGradientNorm<1.0e-10? (ADScalar)0.5 : (massGradient(0)*pressureDirs(0, 0) + massGradient(1)*pressureDirs(1, 0))/massGradientNorm) << std::endl;*/

    //const scalar gamma = 1.0 - (massGradientNorm<1.0e-10? (scalar)0.5 : (massGradient(0)*pressureDirs(0, 0) + massGradient(1)*pressureDirs(1, 0))/massGradientNorm);

    //const scalar velnorm = velocity.norm();
    //Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> Q = massDirs;
    //Q(0, 0) = 1.0; Q(1, 0) = 1.0;
    //Q.col(0) /= Q.col(0).norm();
    //Q(0, 1) = Q(1, 0); Q(1, 1) = -Q(0, 0);

    /*scalar eigscale = 0.5*pressure/energy + 0.5;
    //const scalar para = 1.0-0.5*velocity.dot(velocity)/energy;
    const scalar para = pressure/energy;
    //std::cout << "energy: " << energy << std::endl;
    //std::cout << para << std::endl;
    eigscale = para*0.5 + (1.0-para)*eigscale;*/
    const double eigscale = 0.5;
    Eigen::Matrix<scalar, MacroscaleInformation::dim, 1> eigs;
    eigs(0) = 2.0*eigscale*pressure;
    eigs(1) = 2.0*(1.0-eigscale)*pressure;

    Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> Q = Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>::Identity(MacroscaleInformation::dim, MacroscaleInformation::dim);

    Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = Q*eigs.asDiagonal()*Q.transpose();
    //Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = Q*eigs.asDiagonal()*Q.transpose();
    //Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = pressureDirs*eigs.asDiagonal()*pressureDirs.transpose();
    //Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = massDirs*eigs.asDiagonal()*massDirs.transpose();

    /*for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
      for( std::size_t j=0; j<MacroscaleInformation::dim; ++j ) {
        std::cout << cov_(i, j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;*/

    Eigen::Matrix<scalar, MacroscaleInformation::dim, 1> gamma = Eigen::Matrix<scalar, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);
    //gamma(0) = velocity(0);
    //gamma(1) = velocity(1);
    //gamma(1) = velocity(1)*velocity(1)*massGradient(1)/(std::abs(massGradient.norm())<1.0e-10 ? 1.0 : massGradient.norm());
    //gamma = velJac*massGradient;
    //gamma *= 10.0*mass*energy;
    //gamma(0) = mass*massGradient(0);
    //gamma(1) = mass*massGradient(1);
    //gamma(1) = mass*energy;
    //gamma(0) = mass*pressure*(pressureGradient(1));
    //gamma(1) = mass*pressure*(-pressureGradient(0));
    //gamma(0) = mass*pressure*(massGradNorm*massGradient(0) + pressureGradNorm*pressureGradient(1));
    //gamma(1) = mass*pressure*(massGradNorm*massGradient(1) - pressureGradNorm*pressureGradient(0));
    //gamma(0) = mass*pressureGradNorm*pressureGradient(0);
    //gamma(1) = mass*pressureGradNorm*pressureGradient(1);

    Eigen::Matrix<scalar, MacroscaleInformation::dim, 4> flux;
    flux.col(0) = mass*velocity;
    flux.col(3) = mass*(gamma + energy*velocity + cov*velocity);

    flux.block(0, 1, MacroscaleInformation::dim, MacroscaleInformation::dim) = mass*(velocity*velocity.transpose() + cov);

    /*flux(0, 1) = mass*(velocity(0)*velocity(0) + pressure);
    flux(0, 2) = mass*velocity(0)*velocity(1);

    flux(1, 1) = mass*velocity(1)*velocity(0);
    flux(1, 2) = mass*(velocity(1)*velocity(1) + pressure);*/
    std::cout << "KineticModel Flux Matrix Implementation" << std::endl;
    return flux;
  }

  /**
  @param[in] loc The macro-scale location
  @param[in] mass The macro-scale mass density
  @param[in] velocity The macro-scale (expected) velocity
  @param[in] energy The macro-scale expected energy
  \return The flux matrix for the conservation equations
  */
  virtual Eigen::Matrix<double, MacroscaleInformation::dim, 4> FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const;

  /**
  @param[in] loc The macro-scale location
  @param[in] mass The macro-scale mass density
  @param[in] velocity The macro-scale (expected) velocity
  @param[in] energy The macro-scale expected energy
  \return The flux matrix for the conservation equations
  */
  virtual Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 4> FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, ADScalar const mass, Eigen::Ref<const Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 1> > const& velocity, ADScalar const energy) const;


  /// Get the expected external acceleration at a macro-scale point
  /**
  @param[in] loc The macro-scale location
  @param[in] velocity The expected velocity
  @param[in] time The macro-scale time
  \return The expected external acceleration
  */
  virtual std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double> ExpectedExternalAcceleration(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& loc, Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& velocity, double const time) const;

  /// Get the covariance at a macro-scale point
  /**
  @param[in] loc The macro-scale location
  @param[in] mass The macro-scale mass density
  @param[in] velocity The macro-scale (expected) velocity
  @param[in] energy The macro-scale expected energy
  \return The covariance matrix for the conditional velocity distribution
  */
  virtual Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> Covariance(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const;

  /// Get the covariance at a macro-scale point
  /**
  @param[in] loc The macro-scale location
  @param[in] mass The macro-scale mass density
  @param[in] velocity The macro-scale (expected) velocity
  @param[in] energy The macro-scale expected energy
  \return The skew vector for the conditional velocity distribution
  */
  virtual Eigen::Matrix<double, MacroscaleInformation::dim, 1> Skew(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const;

  /// The initial condition for the mass density
  /**
  Default to return constant \f$1\f$.
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  \return The mass density \f$\mu(\boldsymbol{x})\f$.
  */
  static double InitialMassDensity(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x);

  /// The initial condition for the expected velocity
  /**
  By default, return zero.
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  \return The initial expected velocity \f$\boldsymbol{u}\f$
  */
  static Eigen::Matrix<double, MacroscaleInformation::dim, 1> InitialExpectedVelocity(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x);

  // The initial condition for the expected energy
  /**
  Default to return constant \f$1\f$.
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  \return The expected energy \f$e(\boldsymbol{x})\f$.
  */
  static double InitialExpectedEnergy(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x);

  /// Get the current time
  /**
  The micro-scale models have been updated to this time
  \return The current time
  */
  //double CurrentTime() const;

  /// Get the current time
  /**
  The micro-scale models have been updated to this time
  \return The current time
  */
  //d/ouble& CurrentTime();

  /// Update the micro-scale models
  /**
  By default, the micro-scale models assume the compressible Euler limit, so there are no micro-scale models to update (the conditional density is Guassian). Therefore, this function does nothing.
  @param[in] time The current macro-scale time
  @param[in] macroscaleStates The macro-scale state at grid points
  */
  virtual void UpdateMicroscaleModels(double const time, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& macroscaleStates) const;

  /// Update the fields so that we can compute gradient information at any point
  /**
  @param[in] time The current time
  @param[in] samples The state at locations around the domain
  */
  virtual void UpdateFields(double const time, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples);

  /// The furtherest sqauared distance to a particle method
  const double maxModelDistanceSquared = 0.0;

protected:

  /// The macro-scale states
  std::shared_ptr<spi::Tools::NearestNeighbors> fields;

private:

};

} // namespace KineticEquations
} // namespace spi

#endif
