#ifndef GRAPHLAPLACIAN_HPP_
#define GRAPHLAPLACIAN_HPP_

#include <yaml-cpp/yaml.h>

#include <nanoflann.hpp>

#include <Eigen/Sparse>

#include <MUQ/Modeling/Distributions/RandomVariable.h>

#include <MUQ/SamplingAlgorithms/SampleCollection.h>

#include "spipack/Tools/Kernels/IsotropicKernel.hpp"

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
  such that \f$\mathbb{E}_{\psi}[H] = 0\f$. Let \f$I(i,j)\f$ be the index of the \f$j^{th}\f$ furthest point from the sample \f$\boldsymbol{x}^{(i)}\f$ (note that \f$I(i,0)=i\f$) and let \f$\{ \boldsymbol{x}^{(I(i,j))} \}_{j=0}^{k}\f$ be the \f$k\f$ nearest neighbors. For each sample, define the squared bandwith parameter
  \f{equation*}{
    r_i^2 = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \|^2 \right)}
  \f}

  Define a compact kernel function \f$k: \mathbb{R}^{+} \mapsto \mathbb{R}^{+}\f$ such that \f$k(\theta)\f$ is decreasing and \f$k(\theta) = 0\f$ if \f$\theta \notin [0,1]\f$ (i.e., the support of the kernel is \f$[0,1]\f$). Let
  \f{equation*}{
    k_{l}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(I(i,j))}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \|^2}{2^l r_i r_j} \right)
  \f}



  Define the matrix \f$\boldsymbol{K}\f$ such that the \f$(i^{th}, j^{th})\f$ component is
  \f{equation*}{
    K_{ij} = \frac{k(h^{-2} \|\boldsymbol{x}^{(i)} -\boldsymbol{x}^{(i)}\|^2)}{\sqrt{d^{(i)} d^{(j)}}},
  \f}
  where \f$h\f$ is the bandwidth parameter and \f$d^{(i)} = \sum_{j=1}^{n} k(h^{-2} \|\boldsymbol{x}^{(i)} -\boldsymbol{x}^{(i)}\|^2)\f$. Additionally, define a diagonal matrix \f$\boldsymbol{D}\f$ such that \f$D_{ii} = \sum_{j=1}^{n} K_{ij}\f$. Define the heat matrix and discrete weighted Laplacian
  \f{equation*}{
    \begin{array}{ccc}
      \boldsymbol{P} = \boldsymbol{D}^{-1} \boldsymbol{K} & \mbox{and} & \widetilde{\nabla}_{\psi}^2 = h^{-2} (\boldsymbol{I}-\boldsymbol{P}),
    \end{array}
  \f}
  where \f$\boldsymbol{I}\f$ is the identity matrix.

  <B>Configuration Parameters:</B>
      Parameter Key | Type | Default Value | Description |
      ------------- | ------------- | ------------- | ------------- |
      "NumSamples"   | std::size_t | - | In the case where a random variable is passed to the constructor, we draw \f$n\f$ samples from the distribution.   |
      "MaxLeaf"   | std::size_t | <tt>10</tt> | The maximum leaf size for the kd tree (nanoflann parameter). |
      "KernelOptions"   | YAML::Node | - | The options for the compact kernel function \f$k(\theta)\f$ (spi::Tools::CompactKernel). |
      "Bandwidth"  | double   | <tt>1.0</tt> | The bandwith \f$h\f$ that defines the kernel \f$k_h = k(h^{-1} \theta)\f$.
      "BandwidthIndex"  | int   | <tt>1</tt> | The bandwith index is the maximum value of \f$l\f$ that defines the kernel \f$k_l(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k((2^l r_i r_j)^{-1} \| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2)\f$.
      "EigensolverTol"  | double   | <tt>10^{-5}</tt> | The tolerance for the sparse eigensolver.
      "EigensolverMaxIt"  | std::size_t   | <tt>10^{3}</tt> | The maximum number of iterations for the sparse eigensolver.
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

  /// The squared bandwdwidth \f$h^{2}\f$
  /**
    \return The squared bandwidth \f$h^{2}\f$
  */
  double SquaredBandwidth() const;

  /// Get the \f$i^{th}\f$ point from the point cloud.
  Eigen::Ref<Eigen::VectorXd const> Point(std::size_t const i) const;

  /// Update the kd-tree
  /**
    (re-)Build the kd-tree based on the samples.
  */
  void BuildKDTree();

  /// Find the points that are within radius \f$r\f$ of a given point \f$\boldsymbol{x}\f$
  /**
    Find all of the points \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k} \subseteq \{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$ such that \f$\|\boldsymbol{x}-\boldsymbol{x}^{(j)}\| \leq r\f$, where \f$r>0\f$ is a given radius.
    @param[in] x The given point \f$\boldsymbol{x}\f$
    @param[in] r2 The squared search radius \f$r_i^2 = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \|^2 \right)}\f$
    @param[out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: the squared distance (\f$d^2 = \boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor
  */
  void FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& x, double const r2, std::vector<std::pair<std::size_t, double> >& neighbors) const;

  /// Find the \f$k\f$ points that are closest to a given point \f$\boldsymbol{x}\f$
  /**
    Find the closest \f$k\f$ points \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k} \subseteq \{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$.
    @param[in] x The given point \f$\boldsymbol{x}\f$
    @param[in] k The number of nearest neighbors \f$k\f$
    @param[out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: the squared distance (\f$d^2 = \boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor
    \return The squared distance to the furtherest point from \f$\boldsymbol{x}\f$
  */
  double FindNeighbors(Eigen::Ref<const Eigen::VectorXd> const& x, std::size_t const k, std::vector<std::pair<std::size_t, double> >& neighbors) const;

  /// Evaluate kernel at neighbors
  /**
    Given a point \f$\boldsymbol{x}\f$ and its nearest neighbors \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k}\f$, compute the kernel function \f$k_h(\|\boldsymbol{x} - \boldsymbol{x}^{(j)}\|^2) = k(h^{-2} \|\boldsymbol{x} - \boldsymbol{x}^{(j)}\|^2)\f$.
    @param[in] x The given point \f$\boldsymbol{x}\f$
    @param[in] h2 The squared bandwidth \f$h^2\f$ that defines the kernel.
    @param[in,out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: (input) the squared distance (\f$\boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor, (output) the kernel evaluation \f$k(\boldsymbol{x}, \boldsymbol{x}^{(j)})\f$
    \return The sum of the kernel evaluations \f$d(\boldsymbol{x}) = \sum_{j=1}^{k} k_h(\|\boldsymbol{x}-\boldsymbol{x}^{(j)}\|^2)\f$
  */
  double EvaluateKernel(Eigen::Ref<const Eigen::VectorXd> const& x, double const h2, std::vector<std::pair<std::size_t, double> >& neighbors) const;

  /**
    @param[in] ind The index of the current point
    @param[in] x The given point \f$\boldsymbol{x}^{(i)}\f$
    @param[in] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: the squared distance (\f$\boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor
    @param[in] bandwidth The bandwidth \f$r_i = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \| \right)}\f$
    @param[out] kernelEval The \f$j^{th}\f$ row corresponds to the \f$j^{th}\f$ neighbor and each column stores the kernel evaluation \f$k_l(\boldsymbol{x}, \boldsymbol{x}^{(j)})\f$ for index \f$l\f$.
  */
  void EvaluateKernel(std::size_t const ind, Eigen::Ref<const Eigen::VectorXd> const& x, std::vector<std::pair<std::size_t, double> >& neighbors, Eigen::Ref<Eigen::VectorXd const> const& bandwidth, Eigen::Ref<Eigen::MatrixXd> kernelEval) const;

  /// Construct the heat matrix \f$\boldsymbol{P}\f$
  /**
    Construct the heat matrix \f$\boldsymbol{P}\f$---this function overwrites any existing information in spi::NumericalSolvers::heatMatrix.
  */
  void ConstructHeatMatrix();

  /// Get a reference to the heat matrix
  const Eigen::Ref<const Eigen::SparseMatrix<double> > HeatMatrix() const;

  /// Compute the eigenvalues of the heat matrix \f$\boldsymbol{P}\f$
  /**
    This function does not recompute \f$\boldsymbol{P}\f$, but it computes the largest eigenvalues of the currently stored value.
    @param[in] neig Compute the <tt>neig</tt> largest eigenvalues
    \return The <tt>neig</tt> largest eigenvalues
  */
  Eigen::VectorXd HeatMatrixEigenvalues(const size_t neig) const;

  /// Solve the weighted Poisson problem
  void SolveWeightedPoisson(Eigen::Ref<Eigen::VectorXd> vec);

  /// Get the solver tolerance for the sparse eigenvalue solver
  /**
    \return The solver tolerance for the sparse eigenvalue solver
  */
  double EigensolverTolerance() const;

  /// Write the samples to file
  /**
    @param[in] filename The name of the file where we are writing the data
    @param[in] dataset The name of the dataset where the samples are stored inside the file (defaults to <tt>"/"</tt>)
  */
  void WriteToFile(std::string const& filename, std::string const& dataset = "/") const;

  /// The bandwidth index parameter
  /**
  \return The bandwith index parameter that defines the kernel \f$k_l(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k((2^l r_i r_j)^{-1} \| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2)\f$.
  */
  std::size_t BandwidthIndex() const;

  /// The range of the bandwidth parameter \f$2^{l}\f$
  /**
  We must choose the bandwidth parameter from the range \f$[2^{l_{min}}, 2^{l_{max}}]\f$.
  \return The range \f$[l_{min}, l_{max}]\f$.
  */
  std::pair<double, double> BandwidthRange() const;

  /// Compute the candidate bandwidth parameters
  /**
  \return The possible bandwidth parameters \f$\epsilon = 2^{l_i}\f$.
  */
  Eigen::VectorXd BandwidthParameterCandidates() const;

private:

  /// Create a sample collection by sampling a random variable
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] n Sample the random variable \f$n\f$ times
  */
  static std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> SampleRandomVariable(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, std::size_t const n);

  /// Interpret the particle locations as a point cloud
  struct PointCloud {
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
    Eigen::Ref<Eigen::VectorXd const> Point(std::size_t const i) const;

    /// The dimension of the state
    /**
      The samples \f$x \sim \psi\f$ are in \f$\mathbb{R}^{d}\f$. This function returns the dimension \f$d\f$ (size of the first sample).
    */
    std::size_t StateDim() const;

    /// Samples from the distribution \f$\psi\f$
    std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> samples;
  };

  /// Initialize the graph Laplacian
  /**
    Must be called by every constructor.
    @param[in] options The options for the graph Laplacian
  */
  void Initialize(YAML::Node const& options);

  /// Construct the heat matrix given kernel evaluation at each point's neighbors
  /**
    @param[in] kernelsum The \f$i^{th}\f$ entry is the sum \f$d^{(i)} = \sum_{j=1}^{n} k_h(\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2)\f$
    @param[in] neighbors Each entry is a vector of nearest neighbors for the corresponding point. Each nearest neighbor is the index of the sample and (input) the kernel evaluation \f$k_h(\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2)\f$ (output) the corrected kernel evaluation \f$k_h(\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2) / \sqrt{d^{(i)} d^{(j)}} \f$
    @param[in] numentries A guess of the number of non-zero entries in the matrix (defaults to zero)
  */
  void ConstructHeatMatrix(Eigen::Ref<const Eigen::VectorXd> const& kernelsum, std::vector<std::vector<std::pair<std::size_t, double> > >& neighbors, std::size_t const numentries = 0);

  /// Compute the eigenvalues of a sparse matrix
  /**
    @param[in] neig The number of eigenvalues we wish to compute
    @param[in] mat We are computing the eigenvalues of this matrix
    @param[in] computeLargest True: (default) compute the <tt>neig</tt> largest eigenvalues, False: compute the <tt>neig</tt> smallest eigenvalues.
    \return The eigenvalues
  */
  Eigen::VectorXd ComputeSparseEigenvalues(std::size_t const neig, Eigen::SparseMatrix<double> const& mat, bool const computeLargest = true) const;

  /// Compute the largest eigenvalues of a sparse matrix
  /**
    @param[in] neig The number of eigenvalues we wish to compute
    @param[in] mat We are computing the eigenvalues of this matrix
  */
  Eigen::VectorXd ComputeLargestSparseEigenvalues(std::size_t const neig, Eigen::SparseMatrix<double> const& mat) const;

  /// Compute the smallest eigenvalues of a sparse matrix
  /**
    @param[in] neig The number of eigenvalues we wish to compute
    @param[in] mat We are computing the eigenvalues of this matrix
  */
  Eigen::VectorXd ComputeSmallestSparseEigenvalues(std::size_t const neig, Eigen::SparseMatrix<double> const& mat) const;

  /// The heat matrix \f$\boldsymbol{P}\f$
  Eigen::SparseMatrix<double> heatMatrix;

  /// The point cloud that we will use to approximate the weighted Laplacian (and solve the weighted Poisson equation)
  const PointCloud cloud;

  /// The nanoflann kd-tree type
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud> NanoflannKDTree;

  /// The nanoflann kd-tree
	NanoflannKDTree kdtree;

  /// The kernel function
  /**
    The kernel function \f$k(\theta) = k(h^{-2} \|\boldsymbol{x}_1-\boldsymbol{x}_2\|^2)\f$.
  */
  std::shared_ptr<spi::Tools::IsotropicKernel> kernel;

  /// The maximum and minimum range for the bandwidth index
  /**
  We choose the bandwidth parameter \f$2^{l}\f$; these parameters determine the [min, max] range for \f$l\f$.
  */
  const std::pair<double, double> bandwidthRange;

  /// The bandwidth step parameter
  /**
  We choose the bandwidth parameter from the range \f$[2^{l_{min}}, 2^{l_{max}}]\f$; we discretize this range to \f$\{2^{l_i}\}_{i=1}^{n}\f$. This is the number of points \f$n\f$.
  */
  const std::size_t numBandwidthSteps;

  /// The bandwidth index that defines the kernel \f$k_l(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k((2^l r_i r_j)^{-1} \| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2)\f$.
  const int bandwidthIndex;

  /// The squared bandwidth parameter \f$h^2\f$
  /**
    This parameter defines the kernel \f$k_h(\theta) = k(h^{-2} \theta^2)\f$.
  */
  const double bandwidth2;

  /// The tolerance for the sparse eigensolver.
  const double eigensolverTol;

  /// The maximum number of iterations for the sparse eigensolver.
  const std::size_t eigensolverMaxIt;

  /// The default values for the spi::NumericalSolvers::GraphLaplacian class.
  struct DefaultParameters {
    static double SquaredBandwidth(YAML::Node const& options);

    /// The maximum leaf size (nanoflann parameter) defaults to \f$10\f$
    inline static const std::size_t maxLeaf = 10;

    /// The maximum and minimum range for the bandwidth index defaults to \f$[-1,1]\f$
    /**
    We choose the bandwidth parameter \f$2^{l}\f$; these parameters determine the [min, max] range for \f$l\f$.
    */
    inline static const std::pair<double, double> bandwidthRange = std::pair<double, double>(-1.0, 1.0);

    /// The default number of bandwidth steps is \f$1\f$.
    inline static const std::size_t numBandwidthSteps = 1;

    /// The bandwidth parameter defaults to \f$1\f$
    inline static const double bandwidth = 1.0;

    /// The bandwidth index that defaults to \f$1\f$
    inline static const int bandwidthIndex = 1;

    /// The tolerance for the sparse eigensolver is \f$10^{-5}\f$
    inline static const double eigensolverTol = 1.0e-5;

    /// The maximum number of iterations for the sparse eigensolver \f$10^{3}\f$
    inline static const std::size_t eigensolverMaxIt = 1.0e3;
  };

  /// Store the default values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
