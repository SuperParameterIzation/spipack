#include "spipack/KineticEquations/KineticModels.hpp"

#include "spipack/ConservationEquations/ExaHyPE/Boltzmann/BoltzmannSolver.h"

using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::KineticEquations;

KineticModels::KineticModels() {}

KineticModels::KineticModels(YAML::Node const& options) :
maxModelDistanceSquared(options["MaxModelDistanceSquared"].as<double>(0.0))
{}

double KineticModels::InitialMassDensity(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) {
  //if( x(0)>0.2 & x(0)<0.8 & x(1)>0.2 & x(1)<0.8 ) { return 1.0; }
  //return 0.1;
  return 1.0;
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> KineticModels::InitialExpectedVelocity(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) { return Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1); }

double KineticModels::InitialExpectedEnergy(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) {
  //return ((x-0.5*Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Ones(MacroscaleInformation::dim, 1)).norm()<0.25? 1.0e-1 : 1.0e2);
  return 1.0;
}

Eigen::Matrix<double, MacroscaleInformation::dim, 4> KineticModels::FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const {
  return KineticModels::FluxMatrixImplmentation(loc, mass, velocity, energy);
}

Eigen::Matrix<KineticModels::ADScalar, MacroscaleInformation::dim, 4> KineticModels::FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, ADScalar const mass, Eigen::Ref<const Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 1> > const& velocity, ADScalar const energy) const {
  return KineticModels::FluxMatrixImplmentation(loc, mass, velocity, energy);
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> KineticModels::Skew(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const {
  return Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);
}

Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> KineticModels::Covariance(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const {
  /*// find the nearest neighbors to this location
  assert(nneighbors>1);
  std::vector<std::pair<std::size_t, double> > neighbors;
  const double radius = (1.0+std::log((double)gridSize))/(double)gridSize;
  const double radius2 = radius*radius;
  macroNeighbors->FindNeighbors(loc, radius2, neighbors);
  if( neighbors.size()==0 ) {
    std::cerr << "ERROR: Could not find any neighbors to compute the covariance at location: " << loc.transpose() << std::endl;
    std::cerr << "\tsquared radius: " << radius2 << " num neighbors: " << neighbors.size() << std::endl;
  }
  assert(neighbors.size()>1);

  // fit the regression
  double sm = 0.0;
  Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim> avgcov = Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim>::Zero(MacroscaleInformation::dim, MacroscaleInformation::dim);
  for( std::size_t i=0; i<neighbors.size(); ++i ) {
    assert(neighbors[i].first<models.size());

    const double weight = std::exp(-1.0/(1.0e-10+1.0-neighbors[i].second/radius2));
    assert(!std::isnan(weight)); assert(!std::isinf(weight));
    sm += weight;
    avgcov += weight*models[neighbors[i].first]->Covariance();
  }
  assert(sm>1.0e-10);
  avgcov /= sm;
  assert(!std::isnan(avgcov(0,0)));

  return avgcov;*/
  // compute the pressure
  const double gamma = 0.4; // some parameter that scales the pressure
  const double pressure = gamma*std::max(0.0, energy - 0.5*velocity.dot(velocity));

  return pressure*Eigen::Matrix<double, MacroscaleInformation::dim, MacroscaleInformation::dim>::Identity();
}

std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double> KineticModels::ExpectedExternalAcceleration(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& loc, Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& velocity, double const time) const {
  // the external velocity
  Eigen::Matrix<double, MacroscaleInformation::dim, 1> veldiff = Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1);
  veldiff(0) = 1.0;

  // gyre
  //veldiff(0) = -std::sin(2.0*M_PI*loc(0))*std::cos(2.0*M_PI*loc(1));
  //veldiff(1) = std::cos(2.0*M_PI*loc(0))*std::sin(2.0*M_PI*loc(1));

  // compute the difference between the velocity and the external velocity
  veldiff -= velocity;

  // compute the external acceleration
  return std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double>(veldiff.array().abs()*veldiff.array(), 0.0);
}

void KineticModels::UpdateMicroscaleModels(double const time, std::shared_ptr<SampleCollection> const& macroscaleStates) const {}

std::shared_ptr<KineticModels> KineticModels::ExaHyPEKineticModels() {
  for( const auto& it : exahype::solvers::RegisteredSolvers ) {
    assert(it);
    auto bolt = dynamic_cast<spiEX_Boltzmann::BoltzmannSolver*>(it);
    if( bolt ) { return bolt->kineticModels; }
  }

  return nullptr;
}

void KineticModels::UpdateFields(double const time, std::shared_ptr<SampleCollection> const& samples) {
  // set the options for the graph laplacian
  YAML::Node options;
  options["Stride"] = 1;

  fields = std::make_shared<NearestNeighbors>(samples, options);
  fields->BuildKDTrees();
}
