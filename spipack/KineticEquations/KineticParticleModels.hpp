#ifndef KINETICPARTICLEMODELS_HPP_
#define KINETICPARTICLEMODELS_HPP_

#include <parcer/Communicator.h>


#include "spipack/KineticEquations/MacroscaleInformation.hpp"
#include "spipack/KineticEquations/KineticModels.hpp"

namespace spi {
namespace KineticEquations {

class KineticParticleModels : public KineticModels {
public:

  /// Static construct method for the kinetic models (parallel implementation)
  /**
  @param[in] rv A random variable that we can sample to compute samples from the initial conditional velocity distribution
  @param[in] options Options for the kinetic model
  @param[in] comm The MPI/Parcer communicator
  */
  template<class CONDITIONAL = ConditionalVelocityDistribution, class KINETIC = KineticParticleModels>
  inline static std::shared_ptr<KINETIC> Construct(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node options, std::shared_ptr<parcer::Communicator> const& comm) {
    const std::size_t gridSize = options["GridSize"].as<std::size_t>();

    // the sample collection of macro-scale locations
    auto samples = std::make_shared<muq::SamplingAlgorithms::SampleCollection>();

    // the number of models
    assert(MacroscaleInformation::dim==2);
    const std::size_t numModels = gridSize*gridSize;

    const int rank = comm->GetRank();
    const int mpiSize = comm->GetSize();
    const std::size_t numModelsPerProc = numModels/mpiSize + (rank<numModels%mpiSize);
    const std::size_t start = rank*numModels/mpiSize + std::min(rank, (int)(numModels%mpiSize));
    const std::size_t end = start + numModelsPerProc;

    // create a grid of points where we define a conditional velocity distribution
    std::vector<std::shared_ptr<ConditionalVelocityDistribution> > models(numModelsPerProc);

    const double dx = 1.0/gridSize;
    options["MaxModelDistanceSquared"] = dx*dx/2.0; // max distance is sqrt((dx/2.0)^2+(dx/2.0)^2)
    for( std::size_t n=0; n<numModels; ++n ) {
      // compute the location of the micro-scale model
      Eigen::Matrix<double, MacroscaleInformation::dim, 1> loc;
      loc(0) = (double)(n/gridSize)/(gridSize-1);
      loc(1) = (double)(n%gridSize)/(gridSize-1);
      assert(loc(0)<1.0+1.0e-10); assert(loc(0)>-1.0e-10);
      assert(loc(1)<1.0+1.0e-10); assert(loc(1)>-1.0e-10);
      loc(0) = dx/2.0 + loc(0)*(1.0-dx);
      loc(1) = dx/2.0 + loc(1)*(1.0-dx);
      assert(samples);
      auto samp = std::make_shared<muq::SamplingAlgorithms::SamplingState>(loc);
      assert(samp);
      samples->Add(std::make_shared<muq::SamplingAlgorithms::SamplingState>(loc));

      if( n<start | n>=end ) { continue; }

      // initial mass density
      const double mu = KINETIC::InitialMassDensity(loc);

      // the initial gradient of the mass density
      const Eigen::Matrix<double, MacroscaleInformation::dim, 1> muGrad = KINETIC::InitialMassDensityGradient(loc);

      // initial velocity
      const Eigen::Matrix<double, MacroscaleInformation::dim, 1> vel = KINETIC::InitialExpectedVelocity(loc);

      // initial velocity divergence
      const double velDiv = KINETIC::InitialExpectedVelocityDivergence(loc);

      // initial macro-scale information
      auto macroscale = std::make_shared<MacroscaleInformation>(mu, (mu<1.0e-10 ? Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1).eval() : (muGrad/mu)).eval(), vel, velDiv, Coordinates::MACRO);

      const std::size_t ind = n-start;
      assert(ind<numModelsPerProc);
      models[ind] = std::make_shared<CONDITIONAL>(loc, rv, macroscale, options["ConditionalVelocityDistribution"]);
    }

    return std::make_shared<KINETIC>(samples, models, options, comm);
  }

  /// Static construct method for the kinetic models
  /**
  @param[in] rv A random variable that we can sample to compute samples from the initial conditional velocity distribution
  @param[in] options Options for the kinetic model
  */
  template<class CONDITIONAL = ConditionalVelocityDistribution, class KINETIC = KineticParticleModels>
  inline static std::shared_ptr<KINETIC> Construct(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node options) {
    const std::size_t gridSize = options["GridSize"].as<std::size_t>();

    // the sample collection of macro-scale locations
    auto samples = std::make_shared<muq::SamplingAlgorithms::SampleCollection>();

    // the number of models
    assert(MacroscaleInformation::dim==2);
    const std::size_t numModels = gridSize*gridSize;

    // create a grid of points where we define a conditional velocity distribution
    std::vector<std::shared_ptr<ConditionalVelocityDistribution> > models(numModels);

    const double dx = 1.0/gridSize;
    options["MaxModelDistanceSquared"] = dx*dx/2.0; // max distance is sqrt((dx/2.0)^2+(dx/2.0)^2)
    for( std::size_t n=0; n<numModels; ++n ) {
      // compute the location of the micro-scale model
      Eigen::Matrix<double, MacroscaleInformation::dim, 1> loc;
      loc(0) = (double)(n/gridSize)/(gridSize-1);
      loc(1) = (double)(n%gridSize)/(gridSize-1);
      assert(loc(0)<1.0+1.0e-10); assert(loc(0)>-1.0e-10);
      assert(loc(1)<1.0+1.0e-10); assert(loc(1)>-1.0e-10);
      loc(0) = dx/2.0 + loc(0)*(1.0-dx);
      loc(1) = dx/2.0 + loc(1)*(1.0-dx);
      assert(samples);
      auto samp = std::make_shared<muq::SamplingAlgorithms::SamplingState>(loc);
      assert(samp);
      samples->Add(std::make_shared<muq::SamplingAlgorithms::SamplingState>(loc));

      // initial mass density
      const double mu = KINETIC::InitialMassDensity(loc);

      // the initial gradient of the mass density
      const Eigen::Matrix<double, MacroscaleInformation::dim, 1> muGrad = KINETIC::InitialMassDensityGradient(loc);

      // initial velocity
      const Eigen::Matrix<double, MacroscaleInformation::dim, 1> vel = KINETIC::InitialExpectedVelocity(loc);

      // initial velocity divergence
      const double velDiv = KINETIC::InitialExpectedVelocityDivergence(loc);

      // initial macro-scale information
      auto macroscale = std::make_shared<MacroscaleInformation>(mu, (mu<1.0e-10 ? Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1).eval() : (muGrad/mu)).eval(), vel, velDiv, Coordinates::MACRO);

      models[n] = std::make_shared<CONDITIONAL>(loc, rv, macroscale, options["ConditionalVelocityDistribution"]);
    }

    return std::make_shared<KINETIC>(samples, models, options);
  }

  KineticParticleModels();

  /**
  @param[in] locs The macro-scale location of each conditional model
  @param[in] models The conditional models at key points in the domain
  @param[in] options Options for the kinetic model
  */
  KineticParticleModels(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& locs, std::vector<std::shared_ptr<ConditionalVelocityDistribution> > const& models, YAML::Node const& options);

  /**
  @param[in] locs The macro-scale location of each conditional model
  @param[in] models The conditional models at key points in the domain
  @param[in] options Options for the kinetic model
  @param[in] comm The MPI/Parcer communicator
  */
  KineticParticleModels(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& locs, std::vector<std::shared_ptr<ConditionalVelocityDistribution> > const& models, YAML::Node const& options, std::shared_ptr<parcer::Communicator> const& comm);

  virtual ~KineticParticleModels();

  /// The initial condition for the mass density gradient
  /**
  Default to return zero.
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  \return The gradient of the mass density \f$\nabla_{\boldsymbol{x}} \mu(\boldsymbol{x})\f$.
  */
  static Eigen::Matrix<double, MacroscaleInformation::dim, 1> InitialMassDensityGradient(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x);

  /// The initial condition for the expected velocity divergence
  /**
  By default, return zero.
  @param[in] x The macro-scale location \f$\boldsymbol{x}\f$
  \return The initial expected velocity divergence \f$\nabla_{\boldsymbol{x}} \cdot \boldsymbol{u}\f$
  */
  static double InitialExpectedVelocityDivergence(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x);

  virtual std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double> ExpectedExternalAcceleration(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& loc, Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& velocity, double const time) const override;

  template<typename scalar>
  inline Eigen::Matrix<scalar, MacroscaleInformation::dim, 4> FluxMatrixImplmentation(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, scalar const mass, Eigen::Ref<const Eigen::Matrix<scalar, MacroscaleInformation::dim, 1> > const& velocity, scalar const energy) const {
    // find the nearest small-scale particle model
    std::vector<std::pair<std::size_t, double> > neighbors;
    macroNeighbors->FindNeighbors(loc, (std::size_t)1, neighbors);
    assert(neighbors.size()==1);

    // compute the covariance matrix at this nearest small-scale model, using the macro
    const Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = Covariance<scalar>(neighbors[0].first);
    /*const scalar pressure = std::max(0.0, energy - 0.5*velocity.dot(velocity));
    const Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = pressure*Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>::Identity(MacroscaleInformation::dim, MacroscaleInformation::dim);*/

    const Eigen::Matrix<scalar, MacroscaleInformation::dim, 1> gamma = Skew<scalar>(neighbors[0].first);
    //const Eigen::Matrix<scalar, MacroscaleInformation::dim, 1> gamma = Eigen::Matrix<scalar, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);

    Eigen::Matrix<scalar, MacroscaleInformation::dim, 4> flux;
    flux.col(0) = mass*velocity;
    flux.col(3) = mass*(gamma + energy*velocity + cov*velocity);

    flux.block(0, 1, MacroscaleInformation::dim, MacroscaleInformation::dim) = mass*(velocity*velocity.transpose() + cov);

    return flux;
  }

  /**
  @param[in] loc The macro-scale location
  @param[in] mass The macro-scale mass density
  @param[in] velocity The macro-scale (expected) velocity
  @param[in] energy The macro-scale expected energy
  \return The flux matrix for the conservation equations
  */
  virtual Eigen::Matrix<double, MacroscaleInformation::dim, 4> FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const override;

  /**
  @param[in] loc The macro-scale location
  @param[in] mass The macro-scale mass density
  @param[in] velocity The macro-scale (expected) velocity
  @param[in] energy The macro-scale expected energy
  \return The flux matrix for the conservation equations
  */
  virtual Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 4> FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, ADScalar const mass, Eigen::Ref<const Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 1> > const& velocity, ADScalar const energy) const override;

  /// Get the grid size of the number of kinetic models
  /**
  The kinetic models are organized on a \f$k \times k \f$ grid.
  \return The grid size \f$k\f$.
  */
  std::size_t GridSize() const;

  /// Get the number of micro-scale models
  /**
  \return The number of micro-scale models
  */
  std::size_t NumMicroscaleModels() const;

  /// Get the \f$i^{th}\f$ micro-scale model
  /**
  @param[in] i The index of the model that we want
  \return A pointer to the \f$i^{th}\f$ model
  */
  std::shared_ptr<ConditionalVelocityDistribution> MicroscaleModel(std::size_t const i) const;

  /// Get an iterator to the first model
  /**
  \return An iterator to the first model
  */
  std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::iterator ModelIteratorBegin();

  /// Get a const iterator to the first model
  /**
  \return A const iterator to the first model
  */
  std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::const_iterator ModelIteratorBegin() const;

  /// Get an iterator to the last model
  /**
  \return An iterator to the last model
  */
  std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::iterator ModelIteratorEnd();

  /// Get a const iterator to the last model
  /**
  \return A const iterator to the last model
  */
  std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::const_iterator ModelIteratorEnd() const;

  /// Update the fields so that we can compute gradient information at any point
  /**
  @param[in] time The current time
  @param[in] samples The state at locations around the domain
  */
  virtual void UpdateFields(double const time, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples) override;

private:

  /// Trap any processors that are not rank 0 so that we can update the conditional densities stored on those processors
  void WorkerLoop() const;

  /// Release the the trapped processors
  void ReleaseWorkers() const;

  /// Get the covariance of a specific model
  template<typename scalar>
  Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim> Covariance(std::size_t const index) const {
    if( index<models.size() ) {
      assert(mpiRank==0);
      return models[index]->Covariance().cast<scalar>();
    }

    assert(comm);
    for( std::size_t r=1; r<mpiSize; ++r ) {
      const std::size_t end = (r+1)*gridSize*gridSize/mpiSize + std::min((r+1), (std::size_t)((gridSize*gridSize)%mpiSize));

      if( index<end ) {
        comm->Send(index, r, covarianceTag);
        const Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> cov = comm->Recv<Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> >(r, covarianceTag);
        return cov.cast<scalar>();
      }
    }

    assert(false);
    return Eigen::Matrix<scalar, MacroscaleInformation::dim, MacroscaleInformation::dim>();
  }

  /// Get the skew of a specific model
  template<typename scalar>
  Eigen::Matrix<scalar, MacroscaleInformation::dim, 1> Skew(std::size_t const index) const {
    if( index<models.size() ) {
      assert(mpiRank==0);
      return models[index]->Skew().cast<scalar>();
    }

    assert(comm);
    for( std::size_t r=1; r<mpiSize; ++r ) {
      const std::size_t end = (r+1)*gridSize*gridSize/mpiSize + std::min((r+1), (std::size_t)((gridSize*gridSize)%mpiSize));

      if( index<end ) {
        comm->Send(index, r, skewTag);
        const Eigen::Matrix<double, MacroscaleInformation::dim, 1> skew = comm->Recv<Eigen::Matrix<double, MacroscaleInformation::dim, 1> >(r, skewTag);
        return skew.cast<scalar>();
      }
    }

    assert(false);
    return Eigen::Matrix<scalar, MacroscaleInformation::dim, 1>();
  }

  /// The first (global) ID for the models on this processor
  std::size_t ModelStartIndex() const;

  const int mpiSize = 1;

  const int mpiRank = 0;

  /// Message tags
  const int releaseTag = 2;
  const int covarianceTag = 3;
  const int skewTag = 4;
  const int accelTag = 5;
  const int updateModelsTag = 6;
  const int macroscaleLocationTag = 7;
  const int macroscaleInfoTag = 8;
  const int doneUpdateTag = 9;

  std::shared_ptr<parcer::Communicator> comm;

  /// Update the micro-scale models
  /**
  @param[in] time The current macro-scale time
  */
  void UpdateMicroscaleModels(double const time) const;

  /// Compute the macro-scale information at a specified location given the macro-scale states
  /**
  @param[in] loc The macro-scale location where we need the information
  */
  std::shared_ptr<MacroscaleInformation> ComputeMacroscaleInformation(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& loc) const;

  /// A vector of all the micro-scale models in the domain
  std::vector<std::shared_ptr<ConditionalVelocityDistribution> > models;

  /// The macro-scale points stored in a nearest neighbors object
  std::shared_ptr<spi::Tools::NearestNeighbors> macroNeighbors;

  /// The size of the grid where we define the conditional velocity distributions (the grid will be \f$n \times n\f$)
  const std::size_t gridSize = 1;

  /// The nearest neighbors to used in the regession that allows us to approximate expectations wrt the velocity distribution everywhere
  const std::size_t nneighbors = 10;

  /// The default parameter values for spi::KineticEquations::KineticModels
  struct DefaultParameters {
    /// The default size of the grid where we define the conditional velocity distributions (the grid will be \f$n \times n\f$) is \f$10\f$
    inline static const std::size_t gridSize = 10;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace KineticEquations
} // namespace spi

#endif
