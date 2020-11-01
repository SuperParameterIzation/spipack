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
  This class uses \f$n\f$ samples \f$\{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$ such that \f$\boldsymbol{x}^{(i)} \sim \psi\f$ to approximate the action of this operator and to solve the weighted Poisson equation (given the right hand side \f$R\f$)
  \f{equation*}{
    \nabla_{\psi}^2 H = R
  \f}
  such that \f$\mathbb{E}_{\psi}[H] = 0\f$.
  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "NumSamples"   | std::size_t | - | In the case where a random variable is passed to the constructor, we draw \f$n\f$ samples from the distribution.   |
      "MaxLeaf"   | std::size_t | 10 | The maximum leaf size for the kd tree (nanoflann parameter). |
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
  GraphLaplacian(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options);

  virtual ~GraphLaplacian() = default;

  /// How many samples are in the collection?
  /**
    \return The number of samples in the sample collection (GraphLaplacian::samples)
  */
  std::size_t NumSamples() const;

  /// Get the max leaf size for the kd stree (nanoflann parameter)
  std::size_t KDTreeMaxLeaf() const;

  /// Get the \f$i^{th}\f$ point from the point cloud.
  const Eigen::VectorXd& Point(std::size_t const i) const;

  /// Update the kd-tree
  /**
    (re-)Build the kd-tree based on the samples.
  */
  void BuildKDTree();

  /// Find the points that are within radius \f$r\f$ of a given point \f$\boldsymbol{x}\f$
  /**
    Find all of the points \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k} \subseteq \{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$ such that \f$\|\boldsymbol{x}-\boldsymbol{x}^{(j)}\| \leq r\f$, where \f$r>0\f$ is a given radius.
    @param[in] x The given point \f$\boldsymbol{x}\f$
    @param[in] r The search radius \f$r\f$
    @param[out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: the squared distance (\f$d^2 = \boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor
  */
  void FindNeighbors(Eigen::VectorXd const& x, double const r, std::vector<std::pair<std::size_t, double> >& neighbors) const;

  /// Find the \f$k\f$ points that are closest to a given point \f$\boldsymbol{x}\f$
  /**
    Find the closest \f$k\f$ points \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k} \subseteq \{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$.
    @param[in] x The given point \f$\boldsymbol{x}\f$
    @param[in] k The number of nearest neighbors \f$k\f$
    @param[out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: the squared distance (\f$d^2 = \boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor
    \return The squared distance to the furtherest point from \f$\boldsymbol{x}\f$
  */
  double FindNeighbors(Eigen::VectorXd const& x, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors) const;

  /// Evaluate kernel at neighbors
  /**
    Given a point \f$\boldsymbol{x}\f$ and its nearest neighbors \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k}\f$, compute the kernel function \f$k_h(\|\boldsymbol{x} - \boldsymbol{x}^{(j)}\|^2) = k(h^{-2} \|\boldsymbol{x} - \boldsymbol{x}^{(j)}\|^2)\f$.
    @param[in] x The given point \f$\boldsymbol{x}\f$
    @param[in] h2 The squared bandwidth \f$h^2\f$ that defines the kernel.
    @param[in,out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: (input) the squared distance (\f$d^2 = \boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor, (output) the kernel evaluation \f$k(\boldsymbol{x}, \boldsymbol{x}^{(j)})\f$
  */
  void EvaluateKernel(Eigen::VectorXd const& x, double const h2, std::vector<std::pair<std::size_t, double> >& neighbors) const;

  /// Evaluate the kernel function
  /**
    The kernel is a decreasing function \f$k: \mathbb{R}^{+} \mapsto \mathbb{R}^{+}\f$ such that \f$k(\theta) = 0\f$ if \f$\theta \notin [0,1]\f$.

    The base class implementation is
    \f{equation*}{
      k(\theta) = \begin{cases}
        1 & \mbox{if } \theta \in [0,1] \\
        0 & \mbox{else.}
      \end{cases}
    \f}
    @param[in] theta The input to the kernel \f$\theta\f$
    \return The kernel evaluation \f$k(\theta)\f$
  */
  virtual double Kernel(double const theta) const;

private:

  /// Create a sample collection by sampling a random variable
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] n Sample the random variable \f$n\f$ times
  */
  static std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> SampleRandomVariable(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::size_t const n);

  /// Interpret the particle locations as a point cloud
  class PointCloud {
  public:

    /// Construct the point cloud given samples from the underlying distribution \f$\psi\f$
    /**
      @param[in] samples Samples from the underlying distribution \f$\psi\f$
    */
    PointCloud(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples);

    virtual ~PointCloud() = default;

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

    /// Get the \f$i^{th}\f$ point from the GraphLaplacian::PointCloud::samples.
    const Eigen::VectorXd& Point(std::size_t const i) const;

    /// The dimension of the state
    /**
      The samples \f$x \sim \psi\f$ are in \f$\mathbb{R}^{d}\f$. This function returns the dimension \f$d\f$ (size of the first sample).
    */
    std::size_t StateDim() const;

  private:
    /// Samples from the distribution \f$\psi\f$
    std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samples;
  };

  /// The point cloud that we will use to approximate the weighted Laplacian (and solve the weighted Poisson equation)
  const PointCloud cloud;

  /// The max leaf for the kd-tree
  const std::size_t maxLeaf;

  /// The nanoflann kd-tree type
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud> NanoflannKDTree;

  /// The nanoflann kd-tree
	NanoflannKDTree kdtree;

  /// The default values for the spi::NumericalSolvers::GraphLaplacian class.
  struct DefaultParameters {
    /// The max leaf for the kd-tree defaults to \f$10\f$.
    inline static const std::size_t maxLeaf = 10;
  };

  /// Store the default values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
