#include "spipack/KineticEquations/ConditionalVelocityDistribution.hpp"

#include <MUQ/Utilities/RandomGenerator.h>

#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <sstream>
#include <iomanip>

using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;
using namespace spi::KineticEquations;

ConditionalVelocityDistribution::ConditionalVelocityDistribution(Eigen::VectorXd const& macroLoc, std::shared_ptr<RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
macroLoc(macroLoc),
samples(std::make_shared<NearestNeighbors>(rv, options["NearestNeighbors"])),
kolmogorov(std::make_shared<KolmogorovOperator>(samples, KolmogorovOptions(options["KolmogorovOptions"].as<YAML::Node>(YAML::Node()), samples->StateDim()))),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
varepsilon(options["NondimensionalParameter"].as<double>(defaults.varepsilon)),
alpha(options["ExternalAccelerationRescaling"].as<double>(defaults.alpha)),
theta(options["TimestepParameter"].as<double>(defaults.theta)),
accelerationNoiseScale(options["AccelerationNoiseScale"].as<double>(defaults.accelerationNoiseScale)),
retuneFreq(options["KolmogorovRetuneFrequency"].as<std::size_t>(defaults.retuneFreq)),
numTimesteps(options["NumTimesteps"].as<std::size_t>(defaults.numTimesteps)),
prevMacroInfo(std::make_shared<RescaledMacroscaleInformation>(initMacroInfo, alpha, this->samples->StateDim())),
normalizingConstant(options["InitialNormalizingConstant"].as<double>(defaults.initialNormalizingConstant)),
filename(options["OutputFilename"].as<std::string>(defaults.filename))
{
  assert(theta>-1.0e-10); assert(theta<1.0+1.0e-10);
  assert(alpha>0.0);
}

ConditionalVelocityDistribution::ConditionalVelocityDistribution(Eigen::VectorXd const& macroLoc, std::shared_ptr<SampleCollection> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
macroLoc(macroLoc),
samples(std::make_shared<NearestNeighbors>(samples, options["NearestNeighbors"])),
kolmogorov(std::make_shared<KolmogorovOperator>(this->samples, KolmogorovOptions(options["KolmogorovOptions"].as<YAML::Node>(YAML::Node()), this->samples->StateDim()))),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
varepsilon(options["NondimensionalParameter"].as<double>(defaults.varepsilon)),
alpha(options["ExternalAccelerationRescaling"].as<double>(defaults.alpha)),
theta(options["TimestepParameter"].as<double>(defaults.theta)),
accelerationNoiseScale(options["AccelerationNoiseScale"].as<double>(defaults.accelerationNoiseScale)),
retuneFreq(options["KolmogorovRetuneFrequency"].as<std::size_t>(defaults.retuneFreq)),
numTimesteps(options["NumTimesteps"].as<std::size_t>(defaults.numTimesteps)),
prevMacroInfo(std::make_shared<RescaledMacroscaleInformation>(initMacroInfo, alpha, this->samples->StateDim())),
normalizingConstant(options["InitialNormalizingConstant"].as<double>(defaults.initialNormalizingConstant)),
filename(options["OutputFilename"].as<std::string>(defaults.filename))
{
  assert(theta>-1.0e-10); assert(theta<1.0+1.0e-10);
  assert(alpha>0.0);
}

ConditionalVelocityDistribution::ConditionalVelocityDistribution(Eigen::VectorXd const& macroLoc, std::shared_ptr<NearestNeighbors> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
macroLoc(macroLoc),
samples(samples),
kolmogorov(std::make_shared<KolmogorovOperator>(samples, KolmogorovOptions(options["KolmogorovOptions"].as<YAML::Node>(YAML::Node()), samples->StateDim()))),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
varepsilon(options["NondimensionalParameter"].as<double>(defaults.varepsilon)),
alpha(options["ExternalAccelerationRescaling"].as<double>(defaults.alpha)),
theta(options["TimestepParameter"].as<double>(defaults.theta)),
accelerationNoiseScale(options["AccelerationNoiseScale"].as<double>(defaults.accelerationNoiseScale)),
retuneFreq(options["KolmogorovRetuneFrequency"].as<std::size_t>(defaults.retuneFreq)),
numTimesteps(options["NumTimesteps"].as<std::size_t>(defaults.numTimesteps)),
prevMacroInfo(std::make_shared<RescaledMacroscaleInformation>(initMacroInfo, alpha, this->samples->StateDim())),
normalizingConstant(options["InitialNormalizingConstant"].as<double>(defaults.initialNormalizingConstant)),
filename(options["OutputFilename"].as<std::string>(defaults.filename))
{
  assert(theta>-1.0e-10); assert(theta<1.0+1.0e-10);
  assert(alpha>0.0);
}

YAML::Node ConditionalVelocityDistribution::KolmogorovOptions(YAML::Node options, std::size_t const dim) {
  options["ManifoldDimension"] = (double)dim;
  const YAML::Node& parameter = options["DensityOptions"];
  if( !parameter ) { options["DensityOptions"] = options; }
  options["DensityOptions"]["ManifoldDimension"] = (double)dim;

  return options;
}

std::size_t ConditionalVelocityDistribution::NumSamples() const { return samples->NumSamples(); }

std::size_t ConditionalVelocityDistribution::StateDim() const { return samples->StateDim(); }

double ConditionalVelocityDistribution::CurrentTime() const { return currentTime; }

double ConditionalVelocityDistribution::ExternalAccelerationRescaling() const { return alpha; }

double ConditionalVelocityDistribution::TimestepParameter() const { return theta; }

std::size_t ConditionalVelocityDistribution::NumTimesteps() const { return numTimesteps; }

Eigen::Ref<Eigen::VectorXd const> ConditionalVelocityDistribution::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return samples->Point(i);
}

Eigen::Ref<Eigen::VectorXd> ConditionalVelocityDistribution::Point(std::size_t const i) {
  assert(i<NumSamples());
  return samples->Point(i);
}

void ConditionalVelocityDistribution::Run(double const nextTime, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo, bool const saveInitialConditions) {
  std::cout << "current time: " << currentTime << std::endl;

  // compute the macro-scale time step
  const double macroDelta = nextTime-currentTime;
  assert(macroDelta>0.0);

  // the micro-scale time step
  const double microDelta = 1.0/numTimesteps;
  double microT = -1.0; // the micro-scale time

  // rescale the macro-scale information into the micro-scale coordinates
  auto finalInfo = std::make_shared<RescaledMacroscaleInformation>(finalMacroInfo, alpha, StateDim());
  auto prevInfo = prevMacroInfo;

  std::cout << "macro delta: " << macroDelta << std::endl;
  std::cout << "micro delta: " << microDelta << std::endl;

  // tune the bandwidth parameter
  if( prevInfo->logMassDensityGrad.norm()>1.0e-10 ) {
    kolmogorov->BuildKDTrees();
    kolmogorov->TuneBandwidthParameter(true);
  }

  // write the initial conditions to file
  if( saveInitialConditions & !filename.empty() ) {
    // compute the acceleration at the initial timestep
    ComputeAcceleration(macroDelta, currentTime, prevInfo);

    WriteToFile(currentTime, prevInfo);
  }

  for( std::size_t t=0; t<numTimesteps; ++t ) {
    std::cout << std::endl;

    // the next micro-scale time
    microT += microDelta;
    const double macroTime = currentTime+(t+1)*microDelta*macroDelta;
    std::cout << "next micro time: " << microT << std::endl;

    std::cout << "time step: " << t << " current time: " << macroTime << std::endl;

    // interpolate the macro-scale information on this this timestep
    auto currInfo = InterpolateMacroscaleInformation(microT, finalInfo);

    // update the normalizing constant
    UpdateNormalizingConstant(macroDelta, microDelta, prevInfo, currInfo);

    std::cout << "normalizing constant: " << normalizingConstant << std::endl;

    // collision step
    std::cout << "varepsilon: " << varepsilon << std::endl;
    if( !std::isinf(varepsilon) ) {
      assert(!std::isnan(varepsilon));
      CollisionStep(macroDelta, microDelta, macroTime, currInfo);
    }

    // rebuild the kd tree based on the new particle configuration
    if( currInfo->logMassDensityGrad.norm()>1.0e-10 ) {
      kolmogorov->BuildKDTrees();
      if( (t+1)%retuneFreq==0 & t!=numTimesteps-1 ) { kolmogorov->TuneBandwidthParameter(true); }
    }

    // update the particle velocities
    ConvectionStep(macroDelta, microDelta, macroTime, currInfo);

    // reset the previous macro-scale information
    prevInfo = currInfo;

    // write the current state to file
    if( !filename.empty() ) {
      WriteToFile(macroTime, currInfo);
    }
  }

  // update the current time
  currentTime = nextTime;
  prevMacroInfo = finalInfo;
}

void ConditionalVelocityDistribution::CollisionStep(double const macroDelta, double const microDelta, double const macroTime, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo) {
  const std::size_t n = NumSamples();

  double time = 0.0;
  while( time<microDelta ) {
    // compute the collision probability for each sample
    Eigen::VectorXd collisionProb(n);
    for( std::size_t i=0; i<n; ++i ) { collisionProb(i) = CollisionRateFunction(macroLoc, Point(i), macroTime); }
    collisionProb *= currInfo->massDensity*macroDelta/(varepsilon*std::pow(alpha, (double)StateDim()));

    // compute the timestep size
    const double deltaPrime = std::min(microDelta-time, 1.0/collisionProb.maxCoeff());
    assert(deltaPrime>0.0);
    time += deltaPrime;

    // update the collision probability based on the timestep size
    collisionProb *= deltaPrime;

    // create a list of the particles that will collide
    std::vector<std::size_t> collide;
    collide.reserve(n);
    for( std::size_t i=0; i<n; ++i ) {
      if( RandomGenerator::GetUniform()<collisionProb(i) ) { collide.push_back(i); }
    }

    // randomly shuffle the colliding particles
    std::random_shuffle(std::begin(collide), std::end(collide));

    // collide the particles
    for( int i=0; i<(int)collide.size()-1; i+=2 ) {
      // the indices of the colliding particles
      const std::size_t p1 = collide[i], p2 = collide[i+1];

      // sample and scale a vector from the unit hypersphere
      Eigen::VectorXd w = SampleUnitHypersphere(Point(p1), Point(p2), macroLoc, macroTime);
      if( w.norm()<1.0e-12 ) { continue; }
      assert(std::abs(w.norm()-1.0)<1.0e-10);
      w *= PostCollisionFunction(Point(p1), Point(p2), w, macroLoc, macroTime);

      // update the velocities
      Point(p1) += w; Point(p2) -= w;
    }
  }
}

double ConditionalVelocityDistribution::CollisionRateFunction(Eigen::Ref<const Eigen::VectorXd> const& x, Eigen::Ref<const Eigen::VectorXd> const& v, double const t) { return 1.0; }

Eigen::VectorXd ConditionalVelocityDistribution::SampleUnitHypersphere(Eigen::Ref<const Eigen::VectorXd> const& v, Eigen::Ref<const Eigen::VectorXd> const& vprime, Eigen::Ref<const Eigen::VectorXd> const& x, double const t) const {
  const Eigen::VectorXd diff = v-vprime;
  const double nrm = diff.norm();
  if( nrm<1.0e-12 ) { return diff; }
  return diff/diff.norm();
}

double ConditionalVelocityDistribution::PostCollisionFunction(Eigen::Ref<const Eigen::VectorXd> const& v, Eigen::Ref<const Eigen::VectorXd> const& vprime, Eigen::Ref<const Eigen::VectorXd> const& w, Eigen::Ref<const Eigen::VectorXd> const& x, double const t) const { return -(v-vprime).norm(); }

void ConditionalVelocityDistribution::ConvectionStep(double const macroDelta, double const microDelta, double const macroTime, std::shared_ptr<RescaledMacroscaleInformation> const& currInfo) {
  const std::size_t n = NumSamples();
  ComputeAcceleration(macroDelta, macroTime, currInfo);

  const double delta = macroDelta*microDelta;
  std::cout << "delta: " << delta << " macroDelta: " << macroDelta << " microDelta: " << microDelta << std::endl;
  for( std::size_t i=0; i<n; ++i ) {
    Point(i) += delta*alpha*currInfo->acceleration.row(i);
  }
}

Eigen::VectorXd ConditionalVelocityDistribution::ExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& vel, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const { return Eigen::VectorXd::Zero(StateDim()); }

void ConditionalVelocityDistribution::WeightedPoissonGradient(std::shared_ptr<RescaledMacroscaleInformation> const& currInfo) const {
  const std::size_t n = NumSamples();
  const std::size_t neigs = kolmogorov->NumEigenvalues();

  // compute the eigendecomposition of the kolmogorov operator
  Eigen::VectorXd S(n), Sinv(n), lambda(neigs);
  Eigen::MatrixXd Qhat(n, neigs);
  kolmogorov->ComputeEigendecomposition(S, Sinv, lambda, Qhat);

  // compte the right hand side for the weighted Poisson problem
  Eigen::VectorXd rhs(n);
  for( std::size_t i=0; i<n; ++i ) {
    // the relative velocity: V+Uhat-U(T)=(v/alpha - Uhat) + Uhat - U(T) dotted with the gradient of the log mass density
    rhs(i) = (Point(i)/alpha - currInfo->velocity).dot(currInfo->logMassDensityGrad);
  }
  rhs -= Eigen::VectorXd::Constant(n, rhs.sum()/(double)n);

  // apply the pseudo-inverse to the rhs (we will need the coeffients of the solution)
  Eigen::VectorXd coeff = kolmogorov->PseudoInverse(rhs, S, lambda, Qhat);
  assert(coeff.size()==neigs);

  // compute the gradient of the solution to the weighted Poisson problem
  currInfo->acceleration = kolmogorov->FunctionGradient(coeff, S, Sinv, lambda, Qhat);
}

void ConditionalVelocityDistribution::ComputeAcceleration(double const macroDelta, double const macroTime, std::shared_ptr<RescaledMacroscaleInformation> const& currInfo) {
  const std::size_t dim = StateDim();
  const std::size_t n = NumSamples();

  // is the forcing nonzero---if it is zero, we don't need to compute the augmented acceleration
  if( currInfo->logMassDensityGrad.norm()>1.0e-10 ) {
    WeightedPoissonGradient(currInfo);
  } else {
    currInfo->acceleration = Eigen::MatrixXd::Zero(n, dim);
  }

  const double scale = accelerationNoiseScale*macroDelta;
  assert(scale>-1.0e-10);
  auto gauss = std::make_shared<Gaussian>(dim);

  // add in the external acceleration
  for( std::size_t i=0; i<n; ++i ) {
    Eigen::VectorXd external = ExternalAcceleration(Point(i), macroLoc, macroTime);
    assert(external.size()==dim);
    if( std::abs(scale)>1.0e-12 ) { external += scale*gauss->Sample(); }

    currInfo->acceleration.row(i) = external.transpose()/alpha - currInfo->acceleration.row(i);
  }
}

void ConditionalVelocityDistribution::UpdateNormalizingConstant(double const macroDelta, double const microDelta, std::shared_ptr<const RescaledMacroscaleInformation> const& prevInfo, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo) {
  const double delta = microDelta*macroDelta;
  const double scale = (1.0 + (1.0-theta)*delta*prevInfo->velocityDivergence)/(1.0 - theta*delta*currInfo->velocityDivergence);
  assert(std::abs(scale)>-1.0e-10);

  normalizingConstant *= scale;
}

std::shared_ptr<ConditionalVelocityDistribution::RescaledMacroscaleInformation> ConditionalVelocityDistribution::InterpolateMacroscaleInformation(double const microT, std::shared_ptr<const ConditionalVelocityDistribution::MacroscaleInformation> const& nextMacroInfo) {
  auto macroInfo = std::make_shared<RescaledMacroscaleInformation>();

  // interpolate the macro-scale information
  macroInfo->massDensity = (microT+1.0)*nextMacroInfo->massDensity - microT*prevMacroInfo->massDensity;
  macroInfo->velocity = (microT+1.0)*nextMacroInfo->velocity - microT*prevMacroInfo->velocity;
  macroInfo->velocityDivergence = (microT+1.0)*nextMacroInfo->velocityDivergence - microT*prevMacroInfo->velocityDivergence;
  macroInfo->logMassDensityGrad = (microT+1.0)*nextMacroInfo->logMassDensityGrad - microT*prevMacroInfo->logMassDensityGrad;

  return macroInfo;
}

std::shared_ptr<const ConditionalVelocityDistribution::MacroscaleInformation> ConditionalVelocityDistribution::MacroscaleInfo() const {
  return std::make_shared<MacroscaleInformation>(prevMacroInfo->massDensity/std::pow(alpha, (double)StateDim()), prevMacroInfo->velocity*alpha, prevMacroInfo->velocityDivergence, prevMacroInfo->logMassDensityGrad/alpha);
}

double ConditionalVelocityDistribution::NormalizingConstant() const { return normalizingConstant; }

ConditionalVelocityDistribution::MacroscaleInformation::MacroscaleInformation() {}

ConditionalVelocityDistribution::MacroscaleInformation::MacroscaleInformation(double const massDensity, Eigen::VectorXd const& velocity, double const velocityDivergence, Eigen::VectorXd const& logMassDensityGrad) :
massDensity(massDensity),
velocity(velocity),
velocityDivergence(velocityDivergence),
logMassDensityGrad(logMassDensityGrad)
{}

ConditionalVelocityDistribution::RescaledMacroscaleInformation::RescaledMacroscaleInformation(std::shared_ptr<const MacroscaleInformation> const& macroInfo, double const alpha, std::size_t const dim) :
MacroscaleInformation(macroInfo->massDensity*std::pow(alpha, (double)dim), macroInfo->velocity/alpha, macroInfo->velocityDivergence, macroInfo->logMassDensityGrad*alpha)
{}

ConditionalVelocityDistribution::RescaledMacroscaleInformation::RescaledMacroscaleInformation() :
MacroscaleInformation()
{}

Eigen::MatrixXd ConditionalVelocityDistribution::Covariance(Eigen::Ref<const Eigen::VectorXd> const& mean) const {
  assert(samples);
  return samples->Covariance(mean);
}

double ConditionalVelocityDistribution::ExpectedEnergy(Eigen::Ref<const Eigen::VectorXd> const& expectedVel) const {
  return ExpectedEnergy(0, NumSamples(), expectedVel);
}

Eigen::VectorXd ConditionalVelocityDistribution::ExpectedExternalAcceleration(Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const {
  return ExpectedExternalAcceleration(0, NumSamples(), x, time);
}

double ConditionalVelocityDistribution::ExpectedEnergy(std::size_t const first, std::size_t const last, Eigen::Ref<const Eigen::VectorXd> const& expectedVel) const {
  assert(first<last);
  const std::size_t n = last-first;

  if( n<6 ) {
    double energy = 0.0;
    for( std::size_t i=first; i<last; ++i ) {
      const Eigen::VectorXd diff = Point(i) - expectedVel;
      energy += diff.dot(diff);
    }
    return 0.5*energy/n;
  }

  // compute recursively
  const std::size_t half = n/2;
  const std::size_t middle = first+half;
  return half/(double)n*ExpectedEnergy(first, middle, expectedVel) + (n-half)/(double)n*ExpectedEnergy(middle, last, expectedVel);
}

Eigen::VectorXd ConditionalVelocityDistribution::ExpectedExternalAcceleration(std::size_t const first, std::size_t const last, Eigen::Ref<const Eigen::VectorXd> const& x, double const time) const {
  assert(first<last);
  const std::size_t n = last-first;

  if( n<6 ) {
    Eigen::VectorXd acc = Eigen::VectorXd::Zero(StateDim());
    for( std::size_t i=first; i<last; ++i ) {
      acc += ExternalAcceleration(Point(i), x, time);
    }
    return acc/n;
  }

  // compute recursively
  const std::size_t half = n/2;
  const std::size_t middle = first+half;
  return half/(double)n*ExpectedExternalAcceleration(first, middle, x, time) + (n-half)/(double)n*ExpectedExternalAcceleration(middle, last, x, time);
}

void ConditionalVelocityDistribution::WriteToFile(double const macroTime, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo, std::string const& dataset) const {
  static std::size_t nwrites = 0;
  std::stringstream ss;
  ss << std::setw((std::size_t)log10(numTimesteps)+10) << std::setfill('0') << nwrites++;
  const std::string file = filename+"-"+ss.str()+".h5";

  // output the collection to file
  const std::string dataset_ = (dataset.at(0)=='/'? dataset : "/"+dataset);
  samples->Samples()->WriteToFile(file, dataset_);

  // create an hdf5 file
  auto hdf5file = std::make_shared<HDF5File>(file);
  hdf5file->WriteMatrix(dataset+"/time", Eigen::VectorXd::Constant(1, macroTime).eval());
  hdf5file->WriteMatrix(dataset+"/acceleration", (currInfo->acceleration*alpha).eval());
  hdf5file->WriteMatrix(dataset+"/expected energy", Eigen::VectorXd::Constant(1, ExpectedEnergy(currInfo->velocity*alpha)).eval());
  hdf5file->WriteMatrix(dataset+"/covariance", Covariance(currInfo->velocity*alpha).eval());
  hdf5file->WriteMatrix(dataset+"/expected velocity", (currInfo->velocity*alpha).eval());
  hdf5file->WriteMatrix(dataset+"/sample mean", samples->Mean());
  hdf5file->WriteMatrix(dataset+"/expected acceleration", ExpectedExternalAcceleration(macroLoc, macroTime));
  hdf5file->Close();
}
