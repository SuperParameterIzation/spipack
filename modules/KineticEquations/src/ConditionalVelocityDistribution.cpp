#include "spipack/KineticEquations/ConditionalVelocityDistribution.hpp"

#include <sstream>
#include <iomanip>

using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::NumericalSolvers;
using namespace spi::KineticEquations;

ConditionalVelocityDistribution::ConditionalVelocityDistribution(std::shared_ptr<RandomVariable> const& rv, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
samples(std::make_shared<NearestNeighbors>(rv, options["NearestNeighbors"])),
kolmogorov(std::make_shared<KolmogorovOperator>(samples, KolmogorovOptions(options["KolmogorovOptions"].as<YAML::Node>(YAML::Node()), samples->StateDim()))),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
alpha(options["ExternalAccelerationRescaling"].as<double>(defaults.alpha)),
theta(options["TimestepParameter"].as<double>(defaults.theta)),
numTimesteps(options["NumTimesteps"].as<std::size_t>(defaults.numTimesteps)),
prevMacroInfo(std::make_shared<RescaledMacroscaleInformation>(initMacroInfo, alpha)),
normalizingConstant(options["InitialNormalizingConstant"].as<double>(defaults.initialNormalizingConstant)),
filename(options["OutputFilename"].as<std::string>(defaults.filename))
{
  assert(theta>-1.0e-10); assert(theta<1.0+1.0e-10);
}

ConditionalVelocityDistribution::ConditionalVelocityDistribution(std::shared_ptr<SampleCollection> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
samples(std::make_shared<NearestNeighbors>(samples, options["NearestNeighbors"])),
kolmogorov(std::make_shared<KolmogorovOperator>(this->samples, KolmogorovOptions(options["KolmogorovOptions"].as<YAML::Node>(YAML::Node()), this->samples->StateDim()))),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
alpha(options["ExternalAccelerationRescaling"].as<double>(defaults.alpha)),
theta(options["TimestepParameter"].as<double>(defaults.theta)),
numTimesteps(options["NumTimesteps"].as<std::size_t>(defaults.numTimesteps)),
prevMacroInfo(std::make_shared<RescaledMacroscaleInformation>(initMacroInfo, alpha)),
normalizingConstant(options["InitialNormalizingConstant"].as<double>(defaults.initialNormalizingConstant)),
filename(options["OutputFilename"].as<std::string>(defaults.filename))
{
  assert(theta>-1.0e-10); assert(theta<1.0+1.0e-10);
}

ConditionalVelocityDistribution::ConditionalVelocityDistribution(std::shared_ptr<NearestNeighbors> const& samples, std::shared_ptr<MacroscaleInformation> const& initMacroInfo, YAML::Node const& options) :
samples(samples),
kolmogorov(std::make_shared<KolmogorovOperator>(samples, KolmogorovOptions(options["KolmogorovOptions"].as<YAML::Node>(YAML::Node()), samples->StateDim()))),
currentTime(options["CurrentTime"].as<double>(defaults.currentTime)),
alpha(options["ExternalAccelerationRescaling"].as<double>(defaults.alpha)),
theta(options["TimestepParameter"].as<double>(defaults.theta)),
numTimesteps(options["NumTimesteps"].as<std::size_t>(defaults.numTimesteps)),
prevMacroInfo(std::make_shared<RescaledMacroscaleInformation>(initMacroInfo, alpha)),
normalizingConstant(options["InitialNormalizingConstant"].as<double>(defaults.initialNormalizingConstant)),
filename(options["OutputFilename"].as<std::string>(defaults.filename))
{
  assert(theta>-1.0e-10); assert(theta<1.0+1.0e-10);
}

YAML::Node ConditionalVelocityDistribution::KolmogorovOptions(YAML::Node options, std::size_t const dim) {
  options["ManifoldDimension"] = (double)dim;
  options["DensityOptions"]["ManifoldDimension"] = (double)dim;

  return options;
}

std::size_t ConditionalVelocityDistribution::NumSamples() const { return samples->NumSamples(); }

double ConditionalVelocityDistribution::CurrentTime() const { return currentTime; }

double ConditionalVelocityDistribution::ExternalAccelerationRescaling() const { return alpha; }

double ConditionalVelocityDistribution::TimestepParameter() const { return theta; }

std::size_t ConditionalVelocityDistribution::NumTimesteps() const { return numTimesteps; }

Eigen::Ref<Eigen::VectorXd const> ConditionalVelocityDistribution::Point(std::size_t const i) const {
  assert(i<NumSamples());
  return samples->Point(i);
}

void ConditionalVelocityDistribution::Run(double const nextTime, std::shared_ptr<const MacroscaleInformation> const& finalMacroInfo) {
  std::cout << "current time: " << currentTime << std::endl;

  // compute the macro-scale time step
  const double macroDelta = nextTime-currentTime;
  assert(macroDelta>0.0);

  // the micro-scale time step
  const double microDelta = 1.0/numTimesteps;
  double microT = -1.0; // the micro-scale time

  // rescale the macro-scale information into the micro-scale coordinates
  auto finalInfo = std::make_shared<RescaledMacroscaleInformation>(finalMacroInfo, alpha);
  auto prevInfo = prevMacroInfo;

  std::cout << "macro delta: " << macroDelta << std::endl;
  std::cout << "micro delta: " << microDelta << std::endl;

  // write the initial conditions to file
  std::stringstream ss;
  ss << std::setw((std::size_t)log10(numTimesteps)+1) << std::setfill('0') << 0;
  WriteToFile(filename+"-"+ss.str()+".h5");

  for( std::size_t t=0; t<numTimesteps; ++t ) {
    std::cout << std::endl;

    // the next micro-scale time
    microT += microDelta;
    std::cout << "next micro time: " << microT << std::endl;

    std::cout << "time step: " << t << " current time: " << currentTime+(t+1)*microDelta*macroDelta << std::endl;

    // interpolate the macro-scale information on this this timestep
    auto currInfo = InterpolateMacroscaleInformation(microT, finalInfo);

    // update the normalizing constant
    UpdateNormalizingConstant(macroDelta, microDelta, prevInfo, currInfo);

    std::cout << "normalizing constant: " << normalizingConstant << std::endl;

    // reset the previous macro-scale information
    prevInfo = currInfo;

    // write the current state to file
    std::stringstream ss;
    ss << std::setw((std::size_t)log10(numTimesteps)+1) << std::setfill('0') << t+1;
    WriteToFile(filename+"-"+ss.str()+".h5");
  }

  // update the current time
  currentTime = nextTime;
  prevMacroInfo = finalInfo;
}

void ConditionalVelocityDistribution::UpdateNormalizingConstant(double const macroDelta, double const microDelta, std::shared_ptr<const RescaledMacroscaleInformation> const& prevInfo, std::shared_ptr<const RescaledMacroscaleInformation> const& currInfo) {
  const double delta = microDelta*macroDelta;
  const double scale = (1.0 + (1.0-theta)*delta*prevInfo->velocityDivergence)/(1.0 - theta*delta*currInfo->velocityDivergence);
  assert(std::abs(scale)>-1.0e-10);

  normalizingConstant *= scale;
}

std::shared_ptr<const ConditionalVelocityDistribution::RescaledMacroscaleInformation> ConditionalVelocityDistribution::InterpolateMacroscaleInformation(double const microT, std::shared_ptr<const ConditionalVelocityDistribution::MacroscaleInformation> const& nextMacroInfo) {
  auto macroInfo = std::make_shared<RescaledMacroscaleInformation>();

  // interpolate the macro-scale information
  macroInfo->velocityDivergence = (microT+1.0)*nextMacroInfo->velocityDivergence - microT*prevMacroInfo->velocityDivergence;

  return macroInfo;
}

std::shared_ptr<const ConditionalVelocityDistribution::MacroscaleInformation> ConditionalVelocityDistribution::MacroscaleInfo() const {
  return std::make_shared<MacroscaleInformation>(prevMacroInfo->velocityDivergence);
}

double ConditionalVelocityDistribution::NormalizingConstant() const { return normalizingConstant; }

ConditionalVelocityDistribution::MacroscaleInformation::MacroscaleInformation() {}

ConditionalVelocityDistribution::MacroscaleInformation::MacroscaleInformation(double const velocityDivergence) :
velocityDivergence(velocityDivergence)
{}

ConditionalVelocityDistribution::RescaledMacroscaleInformation::RescaledMacroscaleInformation(std::shared_ptr<const MacroscaleInformation> const& macroInfo, double const alpha) :
MacroscaleInformation(macroInfo->velocityDivergence)
{}

ConditionalVelocityDistribution::RescaledMacroscaleInformation::RescaledMacroscaleInformation() :
MacroscaleInformation()
{}

void ConditionalVelocityDistribution::WriteToFile(std::string const& file, std::string const& dataset) const {
  // output the collection to file
  const std::string dataset_ = (dataset.at(0)=='/'? dataset : "/"+dataset);
  samples->Samples()->WriteToFile(file, dataset_);
}
