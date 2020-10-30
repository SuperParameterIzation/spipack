#ifndef GRAPHLAPLACIAN_HPP_
#define GRAPHLAPLACIAN_HPP_

#include <yaml-cpp/yaml.h>

#include <nanoflann.hpp>

#include <MUQ/Modeling/Distributions/RandomVariable.h>

#include <MUQ/SamplingAlgorithms/SampleCollection.h>

namespace spi {
namespace NumericalSolvers {

/// Estimate a weighted Laplacian operator using a graph Laplacian
/**
  Let \f$\psi\f$ be a probability distribution (with corresponding probability density function). Define the weighted Laplacian operator (given some function \f$H\f$)
  \f{equation*}{
    \nabla_{\psi}^2 H = \psi^{-1} \nabla \cdot (\psi \nabla f).
  \f}
  This class uses \f$n\f$ samples \f$\{x^{(i)}\}_{i=1}^{n}\f$ such that \f$x^{(i)} \sim \psi\f$ to approximate the action of this operator and to solve the weighted Poisson equation (given the right hand side \f$R\f$)
  \f{equation*}{
    \nabla_{\psi}^2 H = R
  \f}
  such that \f$\mathbb{E}_{\psi}[H] = 0\f$.

  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "NumSamples"   | size_t | - | In the case where a random variable is passed to the constructor, we draw \f$n\f$ samples from the distribution.   |
      "MaxLeaf"   | size_t | 10 | The maximum leaf size for the kd tree (nanoflann parameter). |
*/
class GraphLaplacian {
public:

  /// Construct the graph laplacian by sampling a random variable from \f$\psi\f$
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] options Setup options
  */
  GraphLaplacian(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node const& options);

  /// Construct the graph laplacian given samples from the underlying distribution \f$\psi\f$
  /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] options Setup options
  */
  //GraphLaplacian(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options);

  virtual ~GraphLaplacian() = default;

  /// How many samples are in the collection?
  /**
    \return The number of samples in the sample collection (GraphLaplacian::samples)
  */
  size_t NumSamples() const;

  /// Get the max leaf size for the kd stree (nanoflann parameter)
  //size_t KDTreeMaxLeaf() const;

  /// Get the \f$i^{th}\f$ point from the point cloud.
  const Eigen::VectorXd& Point(size_t const i) const;

private:

  /// Create a sample collection by sampling a random variable
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] n Sample the random variable \f$n\f$ times
  */
  static std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> SampleRandomVariable(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, size_t const n);

  /// Interpret the particle locations as a point cloud
  //class PointCloud {
  //public:

    /// Construct the point cloud given samples from the underlying distribution \f$\psi\f$
    /**
      @param[in] samples Samples from the underlying distribution \f$\psi\f$
    */
    //PointCloud(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples);

    //virtual ~PointCloud() = default;

    /// Get the number of samples (from \f$\psi\f$)
    /**
      \return The number of points in the sample collection (GraphLaplacian::PointCloud::samples)
    */
    //size_t kdtree_get_point_count() const;

    /// Get the \f$i^{th}\f$ component of the \f$p^{th}\f$ point
    /**
      @param[in] p We want to access this particle number
      @param[in] i We want this index of the particle location
      \return The \f$i^{th}\f$ component of the \f$p^{th}\f$ point in GraphLaplacian::PointCloud::samples
    */
    //double kdtree_get_pt(size_t const p, size_t const i) const;

    /// Optional bounding-box computation
    /**
      \return Return <tt>false</tt> to default to a standard bounding box computation loop
    */
    /*template<class BBOX>
    inline bool kdtree_get_bbox(BBOX& bb) const { return false; }*/

    /// Get the \f$i^{th}\f$ point from the GraphLaplacian::PointCloud::samples.
    //const Eigen::VectorXd& Point(size_t const i) const;

    /// The dimension of the state
    /**
      The samples \f$x \sim \psi\f$ are in \f$\mathbb{R}^{d}\f$. This function returns the dimension \f$d\f$ (size of the first sample).
    */
    //size_t StateDim() const;

  //private:
    /// Samples from the distribution \f$\psi\f$
    //std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samples;
  //};

  /// The point cloud that we will use to approximate the weighted Laplacian (and solve the weighted Poisson equation)
  //const PointCloud cloud;

  /// The max leaf for the kd-tree
  //const size_t maxLeaf;

  /// The nanoflann kd-tree type
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud> NanoflannKDTree;

  /// The nanoflann kd-tree
	//NanoflannKDTree kdtree;

  /// The default values for the spi::NumericalSolvers::GraphLaplacian class.
  /*struct DefaultParameters {
    /// The max leaf for the kd-tree defaults to \f$10\f$.
    inline static const size_t maxLeaf = 10;
  };

  /// Store the default values
  inline static const DefaultParameters defaults;*/
};

} // ParticleMethods
} // namespace spi

#endif
