#include "spipack/KineticEquations/ConditionalVelocityDistribution.hpp"

#include <MUQ/Utilities/RandomGenerator.h>

#include <MUQ/Modeling/ODE.h>
#include <MUQ/Modeling/Distributions/Gaussian.h>

#include <sstream>
#include <iomanip>

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::KineticEquations;

ConditionalVelocityDistribution::ConditionalVelocityDistribution(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& macroLoc, std::shared_ptr<RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
macroLoc(macroLoc),
samples(std::make_shared<NearestNeighbors>(rv, options["NearestNeighbors"])),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
varepsilon(options["NondimensionalParameter"].as<double>(defaults.varepsilon)),
prevMacroInfo(initMacroInfo)
{
  ShiftSamples(prevMacroInfo);
  InitializeDirectionalDerivative();
  InitializeKolmogorovOperator(options);

  ExpectedEnergy(initMacroInfo->velocity);
}

ConditionalVelocityDistribution::ConditionalVelocityDistribution(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& macroLoc, std::shared_ptr<SampleCollection> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
macroLoc(macroLoc),
samples(std::make_shared<NearestNeighbors>(samples, options["NearestNeighbors"])),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
varepsilon(options["NondimensionalParameter"].as<double>(defaults.varepsilon)),
prevMacroInfo(initMacroInfo)
{
  ShiftSamples(prevMacroInfo);
  InitializeDirectionalDerivative();
  InitializeKolmogorovOperator(options);

  ExpectedEnergy(initMacroInfo->velocity);
}

ConditionalVelocityDistribution::ConditionalVelocityDistribution(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& macroLoc, std::shared_ptr<NearestNeighbors> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
macroLoc(macroLoc),
samples(samples),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
varepsilon(options["NondimensionalParameter"].as<double>(defaults.varepsilon)),
prevMacroInfo(initMacroInfo)
{
  ShiftSamples(prevMacroInfo);
  InitializeDirectionalDerivative();
  InitializeKolmogorovOperator(options);

  ExpectedEnergy(initMacroInfo->velocity);
}

void ConditionalVelocityDistribution::InitializeDirectionalDerivative() {
  // the number of samples
  const std::size_t n = samples->NumSamples();

  directionalDerivative = Eigen::VectorXd(n);
  for( std::size_t i=0; i<n; ++i ) { directionalDerivative(i) = DirectionalDerivative(samples->Point(i)); }
}

double ConditionalVelocityDistribution::DirectionalDerivative(Eigen::VectorXd const& velocity) const { return 0.0; }

void ConditionalVelocityDistribution::InitializeKolmogorovOperator(YAML::Node const& options) {
  // the number of samples
  const std::size_t n = samples->NumSamples();

  const YAML::Node& kolOptions = options["KolmogorovOptions"].as<YAML::Node>(YAML::Node());
  const YAML::Node& kolParaOpt = kolOptions["BandwidthCostOptimization"].as<YAML::Node>(YAML::Node());

  densityEstimationOptions.put("Stride", std::min(n, (std::size_t)(5*log((double)n))));
  densityEstimationOptions.put("NumNearestNeighbors", kolOptions["NumNearestNeighbors"].as<std::size_t>(10));
  densityEstimationOptions.put("ManifoldDimension", (double)MacroscaleInformation::dim);
  densityEstimationOptions.put("SparsityTolerance", kolOptions["SparsityTolerance"].as<double>(1.0e-4));
  densityEstimationOptions.put("NumThreads", kolOptions["NumThreads"].as<std::size_t>(1));
  densityEstimationOptions.put("BandwidthCostOptimization.SparsityTolerance", kolParaOpt["SparsityTolerance"].as<double>(1.0e-4));
  densityEstimationOptions.put("BandwidthCostOptimization.FTol.AbsoluteTolerance", kolParaOpt["FTol.AbsoluteTolerance"].as<double>(1.0e-2));
  densityEstimationOptions.put("BandwidthCostOptimization.FTol.RelativeTolerance", kolParaOpt["FTol.RelativeTolerance"].as<double>(1.0e-2));
  densityEstimationOptions.put("BandwidthCostOptimization.XTol.AbsoluteTolerance", kolParaOpt["XTol.AbsoluteTolerance"].as<double>(1.0e-2));
  densityEstimationOptions.put("BandwidthCostOptimization.XTol.RelativeTolerance", kolParaOpt["XTol.RelativeTolerance"].as<double>(1.0e-2));
  densityEstimationOptions.put("BandwidthCostOptimization.MaxEvaluations", kolParaOpt["MaxEvaluations"].as<std::size_t>(1000));
  densityEstimationOptions.put("BandwidthCostOptimization.Algorithm", kolParaOpt["Algorithm"].as<std::string>("LBFGS"));
  densityEstimationOptions.put("VariableBandwidth", kolOptions["VariableBandwidth"].as<double>(-0.5));
  densityEstimationOptions.put("NumEigenpairs", kolOptions["NumEigenpairs"].as<std::size_t>(std::min(n, (std::size_t)(5*log((double)n)))));

  kolmogorov = std::make_shared<KolmogorovOperator>(samples->Samples(), densityEstimationOptions);
}

std::size_t ConditionalVelocityDistribution::NumSamples() const { return samples->NumSamples(); }

std::size_t ConditionalVelocityDistribution::StateDim() const { return MacroscaleInformation::dim; }

double ConditionalVelocityDistribution::CurrentTime() const { return currentTime; }

Eigen::Ref<Eigen::Matrix<double, MacroscaleInformation::dim, 1> const> ConditionalVelocityDistribution::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return samples->Point(i);
}

Eigen::Ref<Eigen::Matrix<double, MacroscaleInformation::dim, 1> > ConditionalVelocityDistribution::Point(std::size_t const i) {
  assert(i<NumSamples());
  return samples->Point(i);
}

void ConditionalVelocityDistribution::ShiftSamples(std::shared_ptr<const MacroscaleInformation> const& macroInfo) {
  const Eigen::Matrix<double, MacroscaleInformation::dim, 1> shift = macroInfo->velocity-samples->Mean();
  for( std::size_t i=0; i<NumSamples(); ++i ) { Point(i) += shift; }
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::SampleMean() const { return samples->Mean(); }

void ConditionalVelocityDistribution::Run(double const nextTime, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo) {
  // compute the macro-scale time step
  const double delta = nextTime-currentTime;
  if( delta<1.0e-10 ) { return; } // don't need to update anything

  // the number of samples
  const std::size_t n = samples->NumSamples();

  // store the location of the sample pre-collision
  auto tiltedSamples = std::make_shared<muq::SamplingAlgorithms::SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { tiltedSamples->Add(std::make_shared<SamplingState>(kolmogorov->Point(i))); }

  // the collision step
  CollisionStep(delta, finalMacroInfo);

  directionalDerivative = NystromMethod(directionalDerivative, tiltedSamples);

  //const double startTime = currentTime;
  Eigen::VectorXd density = kolmogorov->EstimateDensity();
  const double ssize = nextTime-currentTime;
  currentTime += ssize;

  // create the samples that move off the x=hat{x} plane
  tiltedSamples = std::make_shared<muq::SamplingAlgorithms::SampleCollection>();
  for( std::size_t i=0; i<n; ++i ) { tiltedSamples->Add(std::make_shared<SamplingState>(kolmogorov->Point(i))); }

  // update the samples (both from the distribution and the tilted plane)
  UpdateSamples(prevMacroInfo, density, ssize, tiltedSamples);

  // shift the samples so the macro scale means match
  ShiftSamples(finalMacroInfo);

  // update the densities
  density = kolmogorov->EstimateDensity();
  auto tiltedDensityPoints = std::make_shared<DensityEstimation>(tiltedSamples, densityEstimationOptions);
  const Eigen::VectorXd tiltedDensity = NystromMethod((1.0+ssize*prevMacroInfo->velocityDiv)*(tiltedDensityPoints->EstimateDensity()), tiltedSamples);

  // update the directional derivative
  directionalDerivative = (tiltedDensity-density)/ssize;
  const double nrm = directionalDerivative.sum();
  if( nrm<1.0e-13 ) { directionalDerivative = Eigen::VectorXd::Zero(n); }
  directionalDerivative *= prevMacroInfo->velocityDiv/nrm;

  // update the current time
  assert(std::abs(currentTime-nextTime)<1.0e-13);
  prevMacroInfo = finalMacroInfo;
  computedExpectedAcceleration = false;
  computedCovariance = false;
  computedSkew = false;
  computedEnergy = false;

  ExpectedEnergy(prevMacroInfo->velocity);

  return;
}

Eigen::VectorXd ConditionalVelocityDistribution::NystromMethod(Eigen::VectorXd const& field, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples) const {
  auto density = std::make_shared<DensityEstimation>(samples, densityEstimationOptions);
  //const Eigen::VectorXd tiltedConditional = (1.0+stepsize*currInfo->velocityDiv)*tiltedDensity->EstimateDensity();

  // rebuild the kd trees since the samples have moved
  density->ResetIndices();
  density->BuildKDTrees();
  kolmogorov->ResetIndices();
  kolmogorov->BuildKDTrees();

  // the number of samples
  const std::size_t n = NumSamples();

  // the number of neighbors used to interpolate between samples
  const std::size_t k = 5;

  std::vector<Eigen::Triplet<double> > entries;
  for( std::size_t i=0; i<n; ++i ) {
    // find the points closest to xi
    std::vector<std::pair<std::size_t, double> > neighbors;
    density->FindNeighbors(kolmogorov->Point(i), k, neighbors);
    assert(neighbors.size()==k);

    for( const auto& it : neighbors ) {
      const double eval = std::exp(-it.second);
      entries.emplace_back(i, it.first, eval);
    }
  }

  // the sum of each row in the kernel matrix
  Eigen::VectorXd rowsum = Eigen::VectorXd::Zero(n);
  for( const auto& entry : entries ) { rowsum(entry.row()) += entry.value(); }

  // create the sparse matrix and interpolate
  Eigen::SparseMatrix<double> kernel(n, n);
  kernel.setFromTriplets(entries.begin(), entries.end());

  return rowsum.array().inverse().matrix().asDiagonal()*kernel*field;
}

void ConditionalVelocityDistribution::UpdateSamples(std::shared_ptr<const MacroscaleInformation> const& currInfo, Eigen::VectorXd const& density, double const stepsize, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& tiltedSamples) {
  // the number of samples
  const std::size_t n = NumSamples();

  // compute the eigendecomposition of the Laplace operator
  Eigen::VectorXd similarity, eigs;
  Eigen::MatrixXd Qhat;
  std::tie(similarity, eigs, Qhat) = kolmogorov->Eigendecomposition(density);

  // compute the right hand side for the kolmogorov operator
  Eigen::VectorXd rhs(n);
  assert(rhs.size()==n);
  for( std::size_t i=0; i<n; ++i ) { rhs(i) = (kolmogorov->Point(i)-currInfo->velocity).dot(currInfo->logMassDensityGrad); }

  // effect acceleration of the tilted samples
  Eigen::VectorXd soln = kolmogorov->KolmogorovProblemSolution(similarity, eigs, Qhat, rhs);
  Eigen::MatrixXd accelerationTilted = kolmogorov->GradientVectorField(similarity, eigs, Qhat, soln);

  // update the rhs
  assert(directionalDerivative.size()==n);
  for( std::size_t i=0; i<n; ++i ) { rhs(i) += currInfo->velocityDiv - directionalDerivative(i)/density(i); }

  // effect acceleration of the samples
  soln = kolmogorov->KolmogorovProblemSolution(similarity, eigs, Qhat, rhs);
  Eigen::MatrixXd acceleration = kolmogorov->GradientVectorField(similarity, eigs, Qhat, soln);

  // update using the prescribed acceleration
  const double diffusionNugget = 1.0e-4;
  for( std::size_t i=0; i<n; ++i ) {
    Eigen::Matrix<double, MacroscaleInformation::dim, 1> acc = ExternalAcceleration(Point(i), currentTime);
    accelerationTilted.row(i) += acc;
    acceleration.row(i) += acc;

    // random gaussian
    Eigen::Matrix<double, 1, MacroscaleInformation::dim> vec(MacroscaleInformation::dim, 1);

    for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { vec(i) = RandomGenerator::GetNormal(); }
    tiltedSamples->at(i)->state[0] += stepsize*accelerationTilted.row(i)  + diffusionNugget*vec;

    for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { vec(i) = RandomGenerator::GetNormal(); }
    Point(i) += stepsize*acceleration.row(i) + diffusionNugget*vec;
  }
}

void ConditionalVelocityDistribution::CollisionStep(const double delta, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo) {
  // need to be in the macro-scale coordinate system
  assert(finalMacroInfo->coordinates==Coordinates::MACRO);
  const std::size_t n = NumSamples();

  double time = 0.0;
  while( time<delta ) {
    // the mass density and velocity at this time
    const double dum = delta-time;
    const double massDensity = (time*finalMacroInfo->massDensity + dum*prevMacroInfo->massDensity)/delta;
    const Eigen::VectorXd vel = (time*finalMacroInfo->velocity + dum*prevMacroInfo->velocity)/delta;

    // compute the collision rate function for each sample
    Eigen::VectorXd collisionProb(n);
    const double macroTime = currentTime + time;
    for( std::size_t i=0; i<n; ++i ) { collisionProb(i) = CollisionRateFunction(vel, macroTime); }

    // compute the timestep size
    const double tau = std::min(delta-time, varepsilon/(massDensity*collisionProb.maxCoeff()));

    // compute the collision probability
    collisionProb *= tau*massDensity/varepsilon;

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

      // get the velocities and convert to macro-scale coordinates
      Eigen::Ref<Eigen::Matrix<double, MacroscaleInformation::dim, 1> > v1 = Point(p1);
      Eigen::Ref<Eigen::Matrix<double, MacroscaleInformation::dim, 1> > v2 = Point(p2);

      // sample and scale a vector from the unit hypersphere
      Eigen::VectorXd w = SampleUnitHypersphere(v1, v2, macroTime);
      assert(std::abs(w.norm()-1.0)<1.0e-10);
      w *= PostCollisionFunction(v1, v2, w, macroTime);

      // update the velocities
      v1 += w; v2 -= w;
    }

    time += tau;
  }
}

double ConditionalVelocityDistribution::CollisionRateFunction(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& vi, double const t) const { return 1.0; }

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::SampleUnitHypersphere(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& v, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vprime, double const t) const {
  // random gaussian
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> vec(MacroscaleInformation::dim, 1);
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) { vec(i) = RandomGenerator::GetNormal(); }

  // return the normalized vector
  return vec/vec.norm();
}

double ConditionalVelocityDistribution::PostCollisionFunction(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& v, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vprime, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& w, double const t) const {
  const double we = -w.dot(v-vprime);
  const double e = 0.5*(v.norm()+vprime.norm());
  //const double scale = 1.0e-2;
  //const double gamma = std::exp(scale*prevMacroInfo->massDensity)/(1.0+prevMacroInfo->massDensity*expectedEnergy);
  //std::cout << "gamma: " << gamma << " expected energy: " <<  expectedEnergy << std::endl;
  const double gamma = 1.0;

  return 0.5*(we + std::copysign(1.0, we)*std::sqrt(std::max(0.0, we*we-4.0*(1.0-gamma)*e)));
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::ExternalAcceleration(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& vel, double const time) const {
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> externalVel = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);
  externalVel(0) = 1.0;
  //externalVel(0) = (macroLoc(0)<0.5? -1.0 : 1.0);

  //externalVel -= vel;
  //return externalVel.array().abs()*externalVel.array();
  return externalVel;
}

std::shared_ptr<const MacroscaleInformation> ConditionalVelocityDistribution::MacroscaleInfo() const { return prevMacroInfo; }

Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> ConditionalVelocityDistribution::Covariance() {
  if( !computedCovariance ) {
    assert(samples);
    covariance = samples->Covariance(prevMacroInfo->velocity);
    std::cout << std::endl << std::endl;
    std::cout << covariance << std::endl;
    std::cout << std::endl << std::endl;
    computedCovariance = true;
  }

  return covariance;
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::Skew() {
  if( !computedSkew ) {
    assert(samples);
    skew = Skew(0, samples->NumSamples());
    computedSkew = true;
  }

  return skew;
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::Skew(std::size_t const first, std::size_t const last) const {
  assert(first<=last);
  const std::size_t length = last-first;
  if( length<5 ) {
    Eigen::Matrix<double, MacroscaleInformation::dim, 1> skew = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);
    for( std::size_t i=first; i<last; ++i ) {
      const Eigen::Matrix<double, MacroscaleInformation::dim, 1> diff = samples->Point(i) - prevMacroInfo->velocity;

      skew += diff*diff.dot(diff);
    }
    return skew/(double)length;
  }

  const std::size_t middle = first+length/2;
  const double w1 = (middle-first)/(double)length;
  const double w2 = (last-middle)/(double)length;

  return w1*Skew(first, middle) + w2*Skew(middle, last);
}

double ConditionalVelocityDistribution::ExpectedEnergy(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& expectedVel) {
  if( !computedEnergy ) {
    expectedEnergy = ExpectedEnergy(0, NumSamples(), expectedVel);
    computedEnergy = true;
  }

  return expectedEnergy;
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::ExpectedExternalAcceleration(double const time) {
  if( !computedExpectedAcceleration ) {
    expectedAcceleration = ExpectedExternalAcceleration(0, NumSamples(), time);
    computedExpectedAcceleration = true;
  }

  return expectedAcceleration;
}

double ConditionalVelocityDistribution::ExpectedEnergy(std::size_t const first, std::size_t const last, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& expectedVel) const {
  assert(first<last);
  const std::size_t n = last-first;

  if( n<6 ) {
    double energy = 0.0;
    for( std::size_t i=first; i<last; ++i ) {
      const Eigen::Matrix<double, MacroscaleInformation::dim, 1> diff = Point(i) - expectedVel;
      energy += diff.dot(diff);
    }
    return 0.5*energy/n;
  }

  // compute recursively
  const std::size_t half = n/2;
  const std::size_t middle = first+half;
  return half/(double)n*ExpectedEnergy(first, middle, expectedVel) + (n-half)/(double)n*ExpectedEnergy(middle, last, expectedVel);
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::ExpectedExternalAcceleration(std::size_t const first, std::size_t const last, double const time) const {
  assert(first<last);
  const std::size_t n = last-first;

  if( n<6 ) {
    Eigen::Matrix<double, MacroscaleInformation::dim, 1> acc = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);
    for( std::size_t i=first; i<last; ++i ) {
      Eigen::Matrix<double, MacroscaleInformation::dim, 1> a = ExternalAcceleration(Point(i), time);
      acc += a;
    }
    return acc/n;
  }

  // compute recursively
  const std::size_t half = n/2;
  const std::size_t middle = first+half;
  return half/(double)n*ExpectedExternalAcceleration(first, middle, time) + (n-half)/(double)n*ExpectedExternalAcceleration(middle, last, time);
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> ConditionalVelocityDistribution::MacroscaleLocation() const { return macroLoc; }
