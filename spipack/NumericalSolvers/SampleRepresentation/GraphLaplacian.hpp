#ifndef GRAPHLAPLACIAN_HPP_
#define GRAPHLAPLACIAN_HPP_

#include "spipack/Tools/NearestNeighbors.hpp"

#include "spipack/NumericalSolvers/SampleRepresentation/SampleRepresentation.hpp"

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
      "KernelOptions"   | YAML::Node | - | The options for the compact kernel function \f$k(\theta)\f$ (spi::Tools::CompactKernel). |
      "Bandwidth"  | double   | <tt>1.0</tt> | The bandwith \f$h\f$ that defines the kernel \f$k_h = k(h^{-1} \theta)\f$.
      "EigensolverTol"  | double   | <tt>10^{-5}</tt> | The tolerance for the sparse eigensolver.
      "EigensolverMaxIt"  | std::size_t   | <tt>10^{3}</tt> | The maximum number of iterations for the sparse eigensolver.
*/
class GraphLaplacian : public SampleRepresentation {
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

  /// Do we tune the bandwidth parameter?
  bool TuneBandwidthParameter() const;

  /// Build the kd trees to find nearest neighbors
  void BuildKDTrees();

  /// Evaluate kernel at neighbors
  /**
  Given a point \f$\boldsymbol{x}\f$ and its nearest neighbors \f$\{\boldsymbol{x}^{(j)}\}_{j=1}^{k}\f$, compute the kernel function \f$k_h(\|\boldsymbol{x} - \boldsymbol{x}^{(j)}\|^2) = k(h^{-2} \|\boldsymbol{x} - \boldsymbol{x}^{(j)}\|^2)\f$.
  @param[in] x The given point \f$\boldsymbol{x}\f$
  @param[in] h2 The squared bandwidth \f$h^2\f$ that defines the kernel.
  @param[in,out] neighbors A vector of the nearest neighbors. First: the neighbor's index, Second: (input) the squared distance (\f$\boldsymbol{x} \cdot \boldsymbol{x}^{(j)}\f$) between the point \f$\boldsymbol{x}\f$ and the neighbor, (output) the kernel evaluation \f$k(\boldsymbol{x}, \boldsymbol{x}^{(j)})\f$
  \return The sum of the kernel evaluations \f$d(\boldsymbol{x}) = \sum_{j=1}^{k} k_h(\|\boldsymbol{x}-\boldsymbol{x}^{(j)}\|^2)\f$
  */
  double EvaluateKernel(Eigen::Ref<const Eigen::VectorXd> const& x, double const h2, std::vector<std::pair<std::size_t, double> >& neighbors) const;

  /// Compute the density estimation
  Eigen::VectorXd DensityEstimation() const;

  /**
  @param[in] bandwidth The bandwidth \f$r_i = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \| \right)}\f$
  @param[out] optimal First: The optimal bandwidth parameter, Second: The kernel matrix corresponding to the optimal bandwidth parameter
  \return The sigma parameter
  */
  Eigen::Matrix<double, Eigen::Dynamic, 2> TuneKernelBandwidth(Eigen::Ref<Eigen::VectorXd const> const& bandwidth, std::pair<double, Eigen::SparseMatrix<double> >& optimal) const;

  /**
  @param[in] bandwidth The bandwidth \f$r_i = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \| \right)}\f$
  @param[in] neighbors The nearest neighbors to each point---the kernel is truncated after these neighbors
  @param[out] optimal First: The optimal bandwidth parameter, Second: The kernel matrix corresponding to the optimal bandwidth parameter
  \return The sigma parameter
  */
  Eigen::Matrix<double, Eigen::Dynamic, 2> TuneKernelBandwidth(Eigen::Ref<Eigen::VectorXd const> const& bandwidth, std::vector<std::vector<std::pair<std::size_t, double> > > const& neighbors, std::pair<double, Eigen::SparseMatrix<double> >& optimal) const;

  /**
  @param[in] bandwidthPara The bandwidth parameter \f$\epsilon\f$
  @param[in] bandwidth The bandwidth \f$r_i = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \| \right)}\f$
  @param[out] kernmat The kernel matrix such that \f$K_{ij} = k(\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2/(\epsilon r_i r_j))\f$
  \return Each entry is the sum of the corresponding row
  */
  Eigen::VectorXd KernelMatrix(double const bandwidthPara, Eigen::Ref<Eigen::VectorXd const> const& bandwidth, Eigen::SparseMatrix<double>& kernmat) const;

  /**
  @param[in] bandwidthPara The bandwidth parameter \f$\epsilon\f$
  @param[in] bandwidth The bandwidth \f$r_i = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \| \right)}\f$
  @param[in] neighbors The nearest neighbors to each point---the kernel is truncated after these neighbors
  @param[out] kernmat The kernel matrix such that \f$K_{ij} = k(\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2/(\epsilon r_i r_j))\f$
  \return Each entry is the sum of the corresponding row
  */
  Eigen::VectorXd KernelMatrix(double const bandwidthPara, Eigen::Ref<Eigen::VectorXd const> const& bandwidth, std::vector<std::vector<std::pair<std::size_t, double> > > const& neighbors, Eigen::SparseMatrix<double>& kernmat) const;

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

  /// The bandwidth of each sample
  /**
  \return The squared bandwidth \f$r_i^2 = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \|^2 \right)}\f$
  */
  Eigen::VectorXd SquaredBandwidth() const;

  /// The bandwidth of each sample
  /**
  @param[out] neighbors Each entry is a list of the nearestneighbors for point \f$i\f$. First: the index of the neighbor, Second: the squared distance to the neighbor
  \return The bandwidth \f$r_i^2 = \max_{j \in [1,k]}{\left( \|\boldsymbol{x}^{(i)} - \boldsymbol{x}^{(I(i,j))} \|^2 \right)}\f$
  */
  Eigen::VectorXd SquaredBandwidth(std::vector<std::vector<std::pair<std::size_t, double> > >& neighbors) const;

  /// The range of the bandwidth parameter \f$2^{l}\f$
  /**
  We must choose the bandwidth parameter from the range \f$[2^{l_{min}}, 2^{l_{max}}]\f$.
  \return The range \f$[l_{min}, l_{max}]\f$.
  */
  std::pair<double, double> BandwidthRange() const;

  /// The number of steps in the discretized bandwith range
  /**
  \return The number of steps in the discretized bandwith range
  */
  std::size_t NumBandwidthSteps() const;

  /// Compute the candidate bandwidth parameters
  /**
  \return The possible bandwidth parameters \f$\epsilon = 2^{l_i}\f$.
  */
  Eigen::VectorXd BandwidthParameterCandidates() const;

private:

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

  /// The number of nearest neighbors used to compute the bandwidth
  const std::size_t numNearestNeighbors;

  /// The maximum and minimum range for the bandwidth index
  /**
  We choose the bandwidth parameter \f$2^{l}\f$; these parameters determine the [min, max] range for \f$l\f$.
  */
  const std::pair<double, double> bandwidthRange;

  /// True: tune bandwidth parameter when computing kernel matrices, False: use user-prescribed bandwidth parameter
  const bool tuneBandwidthParameter;

  /// The bandwidth step parameter
  /**
  We choose the bandwidth parameter from the range \f$[2^{l_{min}}, 2^{l_{max}}]\f$; we discretize this range to \f$\{2^{l_i}\}_{i=1}^{n}\f$. This is the number of points \f$n\f$.
  */
  const std::size_t numBandwidthSteps;

  /// The tolerance for the sparse eigensolver.
  const double eigensolverTol;

  /// The maximum number of iterations for the sparse eigensolver.
  const std::size_t eigensolverMaxIt;

  /// The default values for the spi::NumericalSolvers::GraphLaplacian class.
  struct DefaultParameters {
    /// The default number of nearest neighbors for the bandwidth computation is \f$10\f$
    inline static const std::size_t numNearestNeighbors = 10;

    /// The maximum and minimum range for the bandwidth index defaults to \f$[-10,10]\f$
    /**
    We choose the bandwidth parameter \f$2^{l}\f$; these parameters determine the [min, max] range for \f$l\f$.
    */
    inline static const std::pair<double, double> bandwidthRange = std::pair<double, double>(-10.0, 10.0);

    /// The default number of bandwidth steps is \f$10\f$.
    inline static const std::size_t numBandwidthSteps = 10;

    /// By default, do not tune the bandwidth parameter
    inline static const bool tuneBandwidthParameter = false;

    /// The tolerance for the sparse eigensolver is \f$10^{-5}\f$
    inline static const double eigensolverTol = 1.0e-5;

    /// The maximum number of iterations for the sparse eigensolver \f$10^{3}\f$
    inline static const std::size_t eigensolverMaxIt = 1.0e3;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
