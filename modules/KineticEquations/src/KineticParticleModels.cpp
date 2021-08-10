#include "spipack/KineticEquations/KineticParticleModels.hpp"

#include <parcer/Eigen.h>

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

KineticParticleModels::KineticParticleModels(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& locs, std::vector<std::shared_ptr<ConditionalVelocityDistribution> > const& models, YAML::Node const& options, std::shared_ptr<parcer::Communicator> const& comm) :
KineticModels(),
gridSize(options["GridSize"].as<std::size_t>(defaults.gridSize)),
models(models),
mpiSize(comm->GetSize()),
mpiRank(comm->GetRank()),
comm(comm)
{
  if( mpiRank==0 ) {
    // create a nearest neighbor object for the macro-scale locations
    YAML::Node macroNeighborsOptions;
    macroNeighborsOptions["Stride"] = gridSize*gridSize;
    macroNeighbors = std::make_shared<spi::Tools::NearestNeighbors>(locs, macroNeighborsOptions);
    macroNeighbors->BuildKDTrees();
  } else {
    WorkerLoop();
  }
}

KineticParticleModels::~KineticParticleModels() {
  if( (bool)comm & mpiRank==0 ) { ReleaseWorkers(); }
}

void KineticParticleModels::WorkerLoop() const {
  assert(comm);

  // set up a release receive so the master can tell us when we're done
  parcer::RecvRequest continueRequest;
  int releaseComplete;
  bool cont = true;
  comm->Irecv(cont, 0, releaseTag, continueRequest);

  parcer::RecvRequest covRequest;
  parcer::RecvRequest skewRequest;
  parcer::RecvRequest accelRequest;
  parcer::RecvRequest updateModelsRequest;

  const std::size_t start = ModelStartIndex();

  while( true ) {
    if( continueRequest.HasCompleted() ) { if( !cont ) { break; } }

    // check to see if we need a covariance matrix
    const bool needCovariance = comm->Iprobe(0, covarianceTag, covRequest);
    if( needCovariance ) {
      const std::size_t modelID = comm->Recv<std::size_t>(0, covarianceTag)-start;
      assert(modelID>=0);
      assert(modelID<models.size());

      comm->Send(models[modelID]->Covariance(), 0, covarianceTag);
    }

    // check to see if we need a skew vector
    const bool needSkew = comm->Iprobe(0, skewTag, skewRequest);
    if( needSkew ) {
      const std::size_t modelID = comm->Recv<std::size_t>(0, skewTag)-start;
      assert(modelID>=0);
      assert(modelID<models.size());

      comm->Send(models[modelID]->Skew(), 0, skewTag);
    }

    // check to see if we need an acceleration vector
    const bool needAccel = comm->Iprobe(0, accelTag, accelRequest);
    if( needAccel ) {
      std::size_t modelID;
      double time;
      std::tie(modelID, time) = comm->Recv<std::pair<std::size_t, double> >(0, accelTag);
      modelID -= start;
      assert(modelID>=0);
      assert(modelID<models.size());

      comm->Send(models[modelID]->ExpectedExternalAcceleration(time), 0, accelTag);
    }

    // check to see if we need to update the models
    const bool needUpdate = comm->Iprobe(0, updateModelsTag, updateModelsRequest);
    if( needUpdate ) {
      const double time = comm->Recv<double>(0, updateModelsTag);
      UpdateMicroscaleModels(time);
    }
  }
}

std::size_t KineticParticleModels::ModelStartIndex() const { return mpiRank*gridSize*gridSize/mpiSize + std::min(mpiRank, (int)((gridSize*gridSize)%mpiSize)); }

void KineticParticleModels::ReleaseWorkers() const {
  assert(comm);
  const bool falseMsg = false;
  for( int i=1; i<mpiSize; ++i ) { comm->Send(falseMsg, i, releaseTag); }
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

  if( neighbors[0].first<models.size() ) {
    assert(mpiRank==0);
    return std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double>(models[neighbors[0].first]->ExpectedExternalAcceleration(time), neighbors[0].second);
  }

  for( std::size_t r=1; r<mpiSize; ++r ) {
    const std::size_t end = (r+1)*gridSize*gridSize/mpiSize + std::min((r+1), (std::size_t)((gridSize*gridSize)%mpiSize));

    if( neighbors[0].first<end ) {
      comm->Send(std::pair<std::size_t, double>(neighbors[0].first, time), r, accelTag);
      const Eigen::Matrix<double, MacroscaleInformation::dim, 1> accel = comm->Recv<Eigen::Matrix<double, MacroscaleInformation::dim, 1> >(r, accelTag);
      return std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double>(accel, neighbors[0].second);
    }
  }

  assert(false);
  return std::pair<Eigen::Matrix<double, MacroscaleInformation::dim, 1>, double>(Eigen::Matrix<double, MacroscaleInformation::dim, 1>(), neighbors[0].second);
}

std::shared_ptr<MacroscaleInformation> KineticParticleModels::ComputeMacroscaleInformation(Eigen::Matrix<double, MacroscaleInformation::dim, 1> const& loc) const {
  assert(fields);

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
  for( std::size_t r=1; r<mpiSize; ++r ) { comm->Send(time, r, updateModelsTag); }
  UpdateMicroscaleModels(time);

  std::vector<bool> cont(mpiSize-1, true);
  while( true ) {
    for( std::size_t r=1; r<mpiSize; ++r ) {
      parcer::RecvRequest doneUpdateRequest;
      const bool doneUpdate = comm->Iprobe(r, doneUpdateTag, doneUpdateRequest);
      if( doneUpdate ) { cont[r-1] = comm->Recv<bool>(r, doneUpdateTag); }

      if( !cont[r-1] ) { continue; }

      parcer::RecvRequest locRequest;
      const bool needMacroscaleInfo = comm->Iprobe(r, macroscaleLocationTag, locRequest);
      if( needMacroscaleInfo ) {
        Eigen::Matrix<double, MacroscaleInformation::dim, 1> loc = comm->Recv<Eigen::Matrix<double, MacroscaleInformation::dim, 1> >(r, macroscaleLocationTag);
        auto macroInfo = ComputeMacroscaleInformation(loc);
        comm->Send(macroInfo, r, macroscaleInfoTag);
      }
    }

    bool go = false;
    for( std::size_t r=0; r<mpiSize-1; ++r ) {
      if( cont[r] ) {
        go = true;
        break;
      }
    }

    if( !go ) { break; }
  }
}

void KineticParticleModels::UpdateMicroscaleModels(double const time) const {
  // loop through the models and update each one
  //#pragma omp parallel for
  for( const auto& model : models ) {
    assert(model);

    std::cout << "updating model on rank: " << mpiRank << std::endl;

    // compute the macro-scale information at this model location
    std::shared_ptr<MacroscaleInformation> macroInfo;
    if( mpiRank==0 ) {
      parcer::RecvRequest locRequest;
      for( std::size_t r=1; r<mpiSize; ++r ) {
        const bool needMacroscaleInfo = comm->Iprobe(r, macroscaleLocationTag, locRequest);
        if( needMacroscaleInfo ) {
          Eigen::Matrix<double, MacroscaleInformation::dim, 1> loc = comm->Recv<Eigen::Matrix<double, MacroscaleInformation::dim, 1> >(r, macroscaleLocationTag);
          macroInfo = ComputeMacroscaleInformation(loc);
          comm->Send(macroInfo, r, macroscaleInfoTag);
        }
      }

      macroInfo = ComputeMacroscaleInformation(model->MacroscaleLocation());
    } else {
      comm->Send(model->MacroscaleLocation(), 0, macroscaleLocationTag);
      macroInfo = comm->Recv<std::shared_ptr<MacroscaleInformation> >(0, macroscaleInfoTag);
    }
    assert(macroInfo);

    // run the micro-scale model
    model->Run(time, macroInfo);
  }

  if( mpiRank!=0 ) {
    const bool falseMsg = false;
    comm->Send(falseMsg, 0, doneUpdateTag);
  }
}
