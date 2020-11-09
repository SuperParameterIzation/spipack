#ifndef NEARESTNEIGHBORS_HPP_
#define NEARESTNEIGHBORS_HPP_

#include <yaml-cpp/yaml.h>

#include <nanoflann.hpp>

#include <MUQ/Modeling/Distributions/RandomVariable.h>

#include <MUQ/SamplingAlgorithms/SampleCollection.h>

namespace spi {
namespace Tools {

/**
<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"NumSamples"   | std::size_t | - | In the case where a random variable is passed to the constructor, we draw \f$n\f$ samples from the distribution.   |
"MaxLeaf"   | std::size_t | <tt>10</tt> | The maximum leaf size for the kd tree (nanoflann parameter). |
"Stride"   | std::size_t | <tt>NumSamples/100</tt> | Build \f$i \in [0, m]\f$ kd trees that ignore the first \f$i d\f$ samples, where \f$d\f$ is the stride parameter. |
*/
class NearestNeighbors {
public:

  /// Construct the nearest neighbor searcher by sampling a random variable from \f$\psi\f$
  /**
  @param[in] rv The random variable that we wish to sample
  @param[in] options Setup options
  */
  NearestNeighbors(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node const& options);

  /// Construct the nearest neighbor searcher given samples from the underlying distribution \f$\psi\f$
  /**
  @param[in] samples Samples from the underlying distribution \f$\psi\f$
  @param[in] options Setup options
  */
  NearestNeighbors(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options);

  virtual ~NearestNeighbors() = default;

  /// Get the \f$i^{th}\f$ sample
  /**
  @param[in] i We want this sample
  \return Get the \f$i^{th}\f$ sample
  */
  Eigen::Ref<Eigen::VectorXd const> Point(std::size_t const i) const;

  /// The total number of samples
  /**
  \return The number of samples
  */
  std::size_t NumSamples() const;

  /// The state dimension
  /**
  \return The state dimension
  */
  std::size_t StateDim() const;

  /// Build the kd-trees
  void BuildKDTrees() const;

  /// Find the nearest neighbors
  /**
  @param[in] point Find the nearnest neighbors to this point
  @param[in] radius2 The squared search radius
  @param[out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: the squared distance (\f$d^2 = \boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor
  @param[in] lag Ignore the first <tt>lag</tt> samples (defaults to \f$0\f$)
  */
  void FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, double const radius2, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag = 0) const;

  /// Find the nearest neighbors
  /**
  @param[in] point Find the nearnest neighbors to this point
  @param[in] k The number of nearest neighbors \f$k\f$
  @param[out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: the squared distance (\f$d^2 = \boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor
  @param[in] lag Ignore the first <tt>lag</tt> samples (defaults to \f$0\f$)
  */
  void FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag = 0) const;

private:

  /// Create a sample collection by sampling a random variable
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] n Sample the random variable \f$n\f$ times
  */
  static std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> SampleRandomVariable(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::size_t const n);

  /// Initialize the nearest neighbors object
  /**
  @param[in] options Setup options
  */
  void Initialize(YAML::Node const& options);

  /// A point cloud used by nanoflann to find the nearest neighbors
  class Cloud {
  public:
    /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] lag Ignore the first <tt>lag</tt> samples
    */
    Cloud(std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> const& samples, std::size_t const lag);

    virtual ~Cloud() = default;

    /// Get the number of samples (from \f$\psi\f$)
    /**
      \return The number of points in the sample collection (GraphLaplacian::PointCloud::samples)
    */
    std::size_t kdtree_get_point_count() const;

    /// Get the \f$i^{th}\f$ component of the \f$p^{th}\f$ point
    /**
      @param[in] p We want to access this particle number
      @param[in] i We want this index of the particle location
      \return The \f$i^{th}\f$ component of the \f$p^{th}\f$ point in GraphLaplacian::PointCloud::samples
    */
    double kdtree_get_pt(std::size_t const p, std::size_t const i) const;

    /// Optional bounding-box computation
    /**
      \return Return <tt>false</tt> to default to a standard bounding box computation loop
    */
    template<class BBOX>
    inline bool kdtree_get_bbox(BBOX& bb) const { return false; }

    /// Ignore the first <tt>lag</tt> samples
    const std::size_t lag;

  private:

    /// Samples from the distribution \f$\psi\f$
    std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> samples;
  };

  /// The point clouds
	std::vector<Cloud> clouds;

  /// The nanoflann kd-tree type
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Cloud>, Cloud> NanoflannKDTree;

  /// The nanoflann kd-trees
	std::vector<std::shared_ptr<NanoflannKDTree> > kdtrees;

  /// Samples from the distribution \f$\psi\f$
  std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> samples;

  /// The default values for the spi::Tools::NearestNeighbors class.
  struct DefaultParameters {
    /// The maximum leaf size (nanoflann parameter) defaults to \f$10\f$
    inline static const std::size_t maxLeaf = 10;

    /// The stride defaults to 100
    inline static const std::size_t stride = 100;
  };

  /// Store the default values
  inline static const DefaultParameters defaults;
};

} // namespace Tools
} // namespace spi

#endif
