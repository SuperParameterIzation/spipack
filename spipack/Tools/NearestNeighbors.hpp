#ifndef NEARESTNEIGHBORS_HPP_
#define NEARESTNEIGHBORS_HPP_

#include <yaml-cpp/yaml.h>

#include <nanoflann.hpp>

#include <MUQ/Modeling/Distributions/RandomVariable.h>

#include <MUQ/SamplingAlgorithms/SampleCollection.h>

namespace spi {
namespace Tools {

/// Find the nearest neighbors for points in sample collection
/**
Given \f$\{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$ samples from the distribution \f$\psi\f$, this class finds the nearest neighbors \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k}\f$ to a given point \f$\boldsymbol{x}\f$.

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"NumSamples"   | <tt>std::size_t</tt> | - | In the case where a random variable is passed to the constructor, we draw \f$n\f$ samples from the distribution.   |
"MaxLeaf"   | <tt>std::size_t</tt> | <tt>10</tt> | The maximum leaf size for the kd tree (nanoflann parameter). |
"Stride"   | <tt>std::size_t</tt> | <tt>NumSamples/100</tt> | Build \f$i \in [0, m]\f$ kd trees that ignore the first \f$i d\f$ samples, where \f$d\f$ is the stride parameter. |
"NumThreads"   | <tt>std::size_t</tt> | <tt>1</tt> | The number of <tt>openMP</tt> threads available to this object. |
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

  /// The number of <tt>openMP</tt> threads this object can use
  /**
  \return The number of <tt>openMP</tt> threads this object can use
  */
  std::size_t NumThreads() const;

  /// The sample collection
  std::shared_ptr<const muq::SamplingAlgorithms::SampleCollection> Samples() const;

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
  \return The smallest radius that contains all of the points
  */
  double FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& point, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors, std::size_t const& lag = 0) const;

  /// Find the squared bandwidth parameter for each sample
  /**
  The squared bandwidth is defined as \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$ for each sample, where \f$I(i,j)\f$ is the index of the \f$j^{th}\f$ closest neighbor---note that \f$I(i,0) = i\f$.
  @param[in] k The number of nearest neighbors (\f$k\f$)
  @param[out] neighbors Each entry cooresponds to the vector nearest neighbors of the \f$i^{th}\f$ sample. Each neighbor is stored as a pair. First: the index \f$I(i,j)\f$ of the neighbor; second: the squared distance \f$\| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  \return A vector such that the \f$i^{th}\f$ entry is the bandwidth \f$r_i\f$
  */
  Eigen::VectorXd SquaredBandwidth(std::size_t const k, std::vector<std::vector<std::pair<std::size_t, double> > >& neighbors) const;

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

  /// The number of <tt>openMP</tt> threads available to this object.
  const std::size_t numThreads;

  /// The default values for the spi::Tools::NearestNeighbors class.
  struct DefaultParameters {
    /// The maximum leaf size (nanoflann parameter) defaults to \f$10\f$
    inline static const std::size_t maxLeaf = 10;

    /// The default number of <tt>openMP</tt> threads is \f$1\f$.
    inline static const std::size_t numThreads = 1;
  };

  /// Store the default values
  inline static const DefaultParameters defaults;
};

} // namespace Tools
} // namespace spi

#endif
