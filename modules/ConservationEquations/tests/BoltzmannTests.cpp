#include <gtest/gtest.h>

#include <MUQ/Utilities/RandomGenerator.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include "spipack/KineticEquations/KineticParticleModels.hpp"

#include "spipack/ConservationEquations/Boltzmann.hpp"

using namespace muq::Modeling;
using namespace spi::KineticEquations;
using namespace spi::ConservationEquations;

class BoltzmannTests : public::testing::Test {
  /// Set up information to test the nearest neighbor construction
  virtual void SetUp() override {}

  /// Make sure everything is constructed correctly
  virtual void TearDown() override {}
};

TEST_F(BoltzmannTests, ConstructCompressibleEulerLimit) {
  YAML::Node options;
  options["SolverMode"] = "CompressibleEulerLimit";

  // construct the kinetic models
  auto kinetic = std::make_shared<KineticModels>();
  assert(kinetic);

  // the kinetic models are not yet registered
  EXPECT_FALSE(KineticModels::ExaHyPEKineticModels());

  // create the solver
  Boltzmann solver(kinetic, options);

  // make sure the kinetic models are registered
  EXPECT_TRUE(KineticModels::ExaHyPEKineticModels());

  // check that we are using compressible Euler
  EXPECT_EQ(solver.GetSolverMode(), Boltzmann::SolverMode::CompressibleEulerLimit);

  solver.Run();
}

TEST_F(BoltzmannTests, BoltzmannParticleMethods) {
  //Boltzmann::Initialize();

  // the number of particles per micro-scale models
  const std::size_t n = 1000;

  // the number of timesteps per micro-scale model run
  const std::size_t numTimesteps = 5;

  // number of micro-scale models
  std::size_t gridSize = 10;

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = 1;

  // options for the bandwidth parameter optimization within the Kolmogorov problem
  YAML::Node kolParaOptimization;
  kolParaOptimization["SparsityTolerance"] = 1.0e-2;
  kolParaOptimization["FTol.AbsoluteTolerance"] = 1.0e-2;
  kolParaOptimization["FTol.RelativeTolerance"] = 1.0e-2;
  kolParaOptimization["XTol.AbsoluteTolerance"] = 1.0e-2;
  kolParaOptimization["XTol.RelativeTolerance"] = 1.0e-2;
  kolParaOptimization["MaxEvaluations"] = 1000;
  kolParaOptimization["Algorithm"] = "LBFGS";

  // options for the Kolmogorov operator
  YAML::Node kolOptions;
  kolOptions["EigensolverTolerance"] = 1.0e-5;
  kolOptions["TruncationTolerance"] = -std::log(1.0e-4);
  kolOptions["NumNearestNeighbors"] = std::min((std::size_t)25, n/2);
  kolOptions["NumEigenvalues"] = (std::size_t)(5*log((double)n));
  kolOptions["BandwidthCostOptimization"] = kolParaOptimization;

  // options for the conditional velocity distribution
  YAML::Node conditionalOptions;
  conditionalOptions["NearestNeighbors"] = nnOptions;
  conditionalOptions["NumTimesteps"] = numTimesteps;
  conditionalOptions["KolmogorovOptions"] = kolOptions;
  conditionalOptions["AccelerationNoiseScale"] = 0.0;
  conditionalOptions["NondimensionalParameter"] = 1.0e-2;

  // set the options for the kinetic models
  YAML::Node options;
  options["GridSize"] = gridSize;
  options["ConditionalVelocityDistribution"] = conditionalOptions;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(MacroscaleInformation::dim)->AsVariable();

  // construct the micro-scale models
  //{
    auto kinetic = KineticParticleModels::Construct(rv, options);
  //}
  //auto kinetic = std::make_shared<KineticModels>();
  assert(kinetic);

  // the kinetic models are not yet registered
  EXPECT_FALSE(KineticModels::ExaHyPEKineticModels());

  // create the solver
  Boltzmann solver(kinetic, options);

  // make sure the kinetic models are registered
  EXPECT_TRUE(KineticModels::ExaHyPEKineticModels());

  // check that we are using compressible Euler
  EXPECT_EQ(solver.GetSolverMode(), Boltzmann::SolverMode::CompressibleEulerLimit);

  solver.Run();
}

/*TEST_F(BoltzmannTests, ConstructParticleMethods) {
  // the number of samples in each micro-scale simluation
  const std::size_t n = 250;

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = 1;

  // create an (ng x ng) grid of points where we store a conditional velocity distribution
  const std::size_t ng = 20;

  YAML::Node conditionalOptions;
  conditionalOptions["NearestNeighbors"] = nnOptions;

  // create the options
  YAML::Node kineticOptions;
  kineticOptions["GridSize"] = ng;
  kineticOptions["ConditionalVelocityDistribution"] = conditionalOptions;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(MacroscaleInformation::dim)->AsVariable();

  // construct the micro-scale models
  auto models = KineticModels::Construct(rv, kineticOptions);
  assert(models);

  YAML::Node options;
  options["SolverMode"] = "ParticleMethods";

  Boltzmann boltzmannEquation(options);
}*/

class BoltzmannTests_CollisionModel : public spi::KineticEquations::ConditionalVelocityDistribution {
public:
  /// Construct the conditional velocity distribution with a linear external forcing
  /**
  @param[in] options Setup options
  */
  inline BoltzmannTests_CollisionModel(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& macroInfo, YAML::Node const& options) : ConditionalVelocityDistribution(loc, rv, macroInfo, options) {}

  virtual ~BoltzmannTests_CollisionModel() = default;

  /**
  \return The options for this simulation
  */
  inline static YAML::Node Options(std::size_t const n, std::size_t const numTimesteps) {
    // options for the nearest neighbor search
    static YAML::Node nnOptions;
    nnOptions["NumSamples"] = n;
    nnOptions["Stride"] = n/5;
    nnOptions["NumThreads"] = 1;

    // options for the bandwidth parameter optimization within the Kolmogorov problem
    static YAML::Node kolParaOptimization;
    kolParaOptimization["SparsityTolerance"] = 1.0e-1;
    kolParaOptimization["FTolAbsoluteTolerance"] = 1.0e-4;

    // options for the Kolmogorov operator
    static YAML::Node kolOptions;
    kolOptions["EigensolverTolerance"] = 1.0e-5;
    kolOptions["SparsityTolerance"] = 1.0e-4;
    kolOptions["NumNearestNeighbors"] = std::min((std::size_t)25, n/2);
    kolOptions["NumEigenvalues"] = (std::size_t)(5*log((double)n));
    kolOptions["BandwidthCostOptimization.SparsityTolerance"] = 1.0e-4;
    kolOptions["BandwidthCostOptimization"] = kolParaOptimization;

    // set the options for the conditional velocity distribution
    static YAML::Node options;
    options["NearestNeighbors"] = nnOptions;
    options["NumTimesteps"] = numTimesteps;
    options["KolmogorovOptions"] = kolOptions;
    options["AccelerationNoiseScale"] = 0.0;
    options["NondimensionalParameter"] = 1.0e-6;

    return options;
  }

  /**
  @param[in] vel The particle velocity \$\boldsymbol{v}\f$
  @param[in] time The macro-scale time
  \return The external acceleration
  */
  virtual Eigen::Matrix<double, MacroscaleInformation::dim, 1> ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vel, double const time) const override {
    Eigen::Matrix<double, MacroscaleInformation::dim, 1> external = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);
    external(0) = 1.0;

    const Eigen::Matrix<double, MacroscaleInformation::dim, 1> diff = external - vel;

    return 5.0*diff.array().abs()*diff.array();
  }

private:
  /*virtual double PostCollisionFunction(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& v, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vprime, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& w, double const t) const override {
    const double nrg = (v.dot(v)+vprime.dot(vprime))/2.0;
    const double we = -w.dot(v-vprime);

    const double gammamin = (4.0*nrg-we*we)/(4.0*nrg);
    const double gamma = std::max(0.75, gammamin);

    return 0.5*(we + std::copysign(1.0, we)*std::sqrt(std::max(0.0, we-4.0*(1.0-gamma)*nrg)));
  }*/
};

/*
TEST(CompressibleEulerLimitTests, RunModel) {
  // the number of samples in each micro-scale simluation
  const std::size_t n = 250;

  // the number of timesteps in each micro-scale simulation
  const std::size_t numTimesteps = 5;

  // create an (ng x ng) grid of points where we store a conditional velocity distribution
  const std::size_t ng = 20;

  // create the options
  YAML::Node options;
  options["GridSize"] = ng;
  options["ConditionalVelocityDistribution"] = BoltzmannTests_CollisionModel::Options(n, numTimesteps);

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(MacroscaleInformation::dim)->AsVariable();

  // construct the micro-scale models
  auto models = KineticParticleModels::Construct<BoltzmannTests_CollisionModel>(rv, options);
  assert(models);

  //YAML::Node options;
  //options["SolverMode"] = "CompressibleEulerLimit";

  //Boltzmann boltzmannEquation(models, options);

  //boltzmannEquation.Run();
}*/
