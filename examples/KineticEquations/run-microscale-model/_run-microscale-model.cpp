#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <MUQ/Utilities/RandomGenerator.h>

#include "spipack/KineticEquations/ConditionalVelocityDistribution.hpp"

class VelocityDistribution : public spi::KineticEquations::ConditionalVelocityDistribution {
public:
  /// Construct the conditional velocity distribution by sampling a random variable from \f$\psi\f$
  /**
  @param[in] macroLoc The location in the macroscale domain
  @param[in] rv The random variable that we wish to sample
  @param[in] initMacroInfo The initial macro-scale information
  @param[in] options Setup options
  */
  inline VelocityDistribution(Eigen::VectorXd const& macroLoc, std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) : ConditionalVelocityDistribution(macroLoc, rv, initMacroInfo, options) {}

  virtual ~VelocityDistribution() = default;

  inline virtual Eigen::VectorXd ExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& vel, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const override {
    Eigen::VectorXd external = Eigen::VectorXd::Zero(StateDim());
    //return external;
    external(0) = 1.0;

    external = external-vel;
    return external.norm()*external;

    //return external.array().abs()*external.array();
  }

private:

  inline virtual Eigen::VectorXd SampleUnitHypersphere(Eigen::Ref<const Eigen::VectorXd> const& v, Eigen::Ref<const Eigen::VectorXd> const& vprime, Eigen::Ref<const Eigen::VectorXd> const& x, double const t) const {
    const double theta = 2.0*M_PI*muq::Utilities::RandomGenerator::GetUniform();
    Eigen::VectorXd w(2);
    w(0) = std::cos(theta); w(1) = std::sin(theta);
    return w;
  }

  inline virtual double PostCollisionFunction(Eigen::Ref<const Eigen::VectorXd> const& v, Eigen::Ref<const Eigen::VectorXd> const& vprime, Eigen::Ref<const Eigen::VectorXd> const& w, Eigen::Ref<const Eigen::VectorXd> const& x, double const t) const override {
    //double gamma = 0.8+0.3*std::sin(2.0*M_PI*t/0.1);
    //double gamma = 0.95;
    double gamma = 1.0;
    const double sigmae = -(v-vprime).dot(w);

    //std::cout << "sige: " << sigmae << std::endl;

    const double e = (v.dot(v)+vprime.dot(vprime))/2.0;
    assert(e>1.0e-10);
    gamma = std::max(gamma, (4.0*e-sigmae*sigmae)/(4.0*e)+1.0e-5*gamma);
    //gamma = std::max(gamma, (4.0*e-sigmae*sigmae)/(4.0*e));
    //std::cout << "gamma: " << gamma << std::endl;

    //std::cout << "sqrt: " << sigmae*sigmae-4.0*(1.0-gamma)*e << std::endl;

    //std::cout << 0.5*sigmae + std::copysign(1.0, sigmae)/2.0*std::sqrt(std::abs(sigmae*sigmae-4.0*(1.0-gamma)*e)) << std::endl;

    // max gamma ensures that the thing inside the square root should be positive--use max(x, 0) to avoid roundoff errors
    return 0.5*sigmae + std::copysign(1.0, sigmae)/2.0*std::sqrt(std::max(0.0, sigmae*sigmae-4.0*(1.0-gamma)*e));
  }
};

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::KineticEquations;

int main(int argc, char **argv) {
  // the dimension of state spaces
  const unsigned int dim = 2;

  // create a standard Gaussian random variable
  auto rv = std::make_shared<Gaussian>(dim)->AsVariable();

  // the number of samples
  const std::size_t n = 1000;

  // the number of timesteps
  const std::size_t numTimesteps = 100;

  // the output filename
  std::string filename = "output/output";

  // the macroscale location
  const Eigen::VectorXd macroLoc = Eigen::VectorXd::Zero(dim);

  // options for the nearest neighbor search
  YAML::Node nnOptions;
  nnOptions["NumSamples"] = n;
  nnOptions["Stride"] = n/5;
  nnOptions["NumThreads"] = omp_get_max_threads();

  // set the options for the conditional velocity distribution
  YAML::Node options;
  options["NearestNeighbors"] = nnOptions;
  options["NumTimesteps"] = numTimesteps;
  options["OutputFilename"] = filename;
  options["NondimensionalParameter"] = 1.0;

  // construct the initial macro-scale information
  const double initMassDensity = 1.0;
  const Eigen::VectorXd initVel = Eigen::VectorXd::Zero(dim);
  const double initVelDiv = 0.0;
  const Eigen::VectorXd initLogMassDensityGrad = Eigen::VectorXd::Zero(dim);
  auto macroInfo = std::make_shared<ConditionalVelocityDistribution::MacroscaleInformation>(initMassDensity, initVel, initVelDiv, initLogMassDensityGrad);

  // create the conditional velocity distribution
  auto distribution = std::make_shared<VelocityDistribution>(macroLoc, rv, macroInfo, options);

  // construct the macro-scale information
  double massDensity = 1.0;
  Eigen::VectorXd vel = Eigen::VectorXd::Zero(dim);
  //vel(0) = 1.0;
  const double velDiv = 0.0;
  Eigen::VectorXd logMassDensityGrad = Eigen::VectorXd::Zero(dim);
  logMassDensityGrad(0) = 0.0;
  auto nextMacroInfo = std::make_shared<ConditionalVelocityDistribution::MacroscaleInformation>(massDensity, vel, velDiv, logMassDensityGrad);

  // run the micro-scale model
  double time = 0.0;
  const double delta = 0.1;
  for( std::size_t i=0; i<5; ++i ) {
    time += delta;
    distribution->Run(time, nextMacroInfo);
  }
}
