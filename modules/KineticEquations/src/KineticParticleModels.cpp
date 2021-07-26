#include "spipack/KineticEquations/KineticParticleModels.hpp"

using namespace muq::SamplingAlgorithms;
using namespace spi::Tools;
using namespace spi::KineticEquations;

KineticParticleModels::KineticParticleModels() : KineticModels() {}

KineticParticleModels::KineticParticleModels(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& locs, std::vector<std::shared_ptr<ConditionalVelocityDistribution> > const& models, YAML::Node const& options) :
KineticModels(),
gridSize(options["GridSize"].as<std::size_t>(defaults.gridSize)),
models(models)
{
  // create a nearest neighbor object for the macro-scale locations
  YAML::Node macroNeighborsOptions;
  macroNeighborsOptions["Stride"] = gridSize*gridSize;
  macroNeighbors = std::make_shared<spi::Tools::NearestNeighbors>(locs, macroNeighborsOptions);
  macroNeighbors->BuildKDTrees();
}

Eigen::Matrix<double, MacroscaleInformation::dim, 1> KineticParticleModels::InitialMassDensityGradient(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) { return Eigen::Matrix<double, MacroscaleInformation::dim, 1>::Zero(MacroscaleInformation::dim, 1); }

double KineticParticleModels::InitialExpectedVelocityDivergence(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& x) { return 0.0; }

std::size_t KineticParticleModels::GridSize() const { return gridSize; }

std::size_t KineticParticleModels::NumMicroscaleModels() const { return models.size(); }

std::shared_ptr<ConditionalVelocityDistribution> KineticParticleModels::MicroscaleModel(std::size_t const i) const { return models[i]; }

std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::iterator KineticParticleModels::ModelIteratorBegin() { return models.begin(); }

std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::const_iterator KineticParticleModels::ModelIteratorBegin() const { return models.begin(); }

std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::iterator KineticParticleModels::ModelIteratorEnd() { return models.end(); }

std::vector<std::shared_ptr<ConditionalVelocityDistribution> >::const_iterator KineticParticleModels::ModelIteratorEnd() const { return models.end(); }

Eigen::Matrix<double, MacroscaleInformation::dim, 4> KineticParticleModels::FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, double const mass, Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& velocity, double const energy) const {
  return KineticParticleModels::FluxMatrixImplmentation(loc, mass, velocity, energy);
}

Eigen::Matrix<KineticModels::ADScalar, MacroscaleInformation::dim, 4> KineticParticleModels::FluxMatrix(Eigen::Ref<const Eigen::Matrix<double, MacroscaleInformation::dim, 1> > const& loc, ADScalar const mass, Eigen::Ref<const Eigen::Matrix<ADScalar, MacroscaleInformation::dim, 1> > const& velocity, ADScalar const energy) const {
  return KineticParticleModels::FluxMatrixImplmentation(loc, mass, velocity, energy);
}

std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double> KineticParticleModels::ExpectedExternalAcceleration(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& loc, Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& velocity, double const time) const {
  // find the nearest small-scale particle model
  std::vector<std::pair<std::size_t, double> > neighbors;
  macroNeighbors->FindNeighbors(loc, (std::size_t)1, neighbors);
  assert(neighbors.size()==1);

  return std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double>(models[neighbors[0].first]->ExpectedExternalAcceleration(time), neighbors[0].second);
}

std::shared_ptr<MacroscaleInformation> KineticParticleModels::ComputeMacroscaleInformation(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& loc) const {
  constexpr std::size_t nNeighs = 5;
  std::vector<std::pair<std::size_t, double> > neighbors;
  fields->FindNeighbors(loc, nNeighs, neighbors);
  assert(neighbors.size()==nNeighs);

  // the Vandermonde matrix and data for linear regression
  Eigen::Matrix<double, nNeighs, MacroscaleInformation::dim+1> vand;
  Eigen::Matrix<double, nNeighs, MacroscaleInformation::dim+2> data;

  // loop through the neighbors
  for( std::size_t i=0; i<nNeighs; ++i ) {
    auto sample =  fields->GetSamplingState(neighbors[i].first);

    // the Vandermonde matrix
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

  // store the macro-scale information
  auto macroInfo = std::make_shared<MacroscaleInformation>();
  macroInfo->coordinates = Coordinates::MACRO;
  macroInfo->massDensity = coeff(0, 0) + coeff.col(0).tail(MacroscaleInformation::dim).dot(loc);
  macroInfo->logMassDensityGrad = coeff.col(0).tail(MacroscaleInformation::dim)/std::max(1.0e-10, macroInfo->massDensity);
  macroInfo->velocityDiv = 0.0;
  for( std::size_t i=0; i<MacroscaleInformation::dim; ++i ) {
    macroInfo->velocity(i) = coeff(0, i+1) + coeff.col(i+1).tail(MacroscaleInformation::dim).dot(loc);
    macroInfo->velocityDiv += coeff(i+1, i+1);
  }

  return macroInfo;
}

void KineticParticleModels::UpdateFields(double const time, std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples) {
  KineticModels::UpdateFields(time, samples);

  // loop through the models
  UpdateMicroscaleModels(time);
}

void KineticParticleModels::UpdateMicroscaleModels(double const time) const {
  // loop through the models and update each one
  #pragma omp parallel for
  for( const auto& model : models ) {
    assert(model);

    // compute the macro-scale information at this model location
    auto macroInfo = ComputeMacroscaleInformation(model->MacroscaleLocation());
    assert(macroInfo);

    // run the micro-scale model
    model->Run(time, macroInfo);

    // rebuild the estimates of the expecations
    //model->ExpectedExternalAcceleration(time);
    //model->Covariance(macroInfo->velocity);
  }
}
