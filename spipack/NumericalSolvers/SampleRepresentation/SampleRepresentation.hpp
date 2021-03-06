#ifndef SAMPLEREPRESENTATION_HPP_
#define SAMPLEREPRESENTATION_HPP_

#include <boost/property_tree/ptree.hpp>

#include <Eigen/Sparse>

#include "spipack/Tools/NearestNeighbors.hpp"

#include "spipack/Tools/Kernels/IsotropicKernel.hpp"

namespace spi {
namespace NumericalSolvers {

/// Represent a distribution \f$\psi\f$ using samples
/**
Let \f$\psi\f$ be a probability distribution, which we represent here using \f$n\f$ samples \f$\{\boldsymbol{x}^{(i)}\}_{i=1}^{n}\f$ such that \f$\boldsymbol{x}^{(i)} \sim \psi\f$.

Let \f$r_i\f$ be a bandwidth  parmeter associated with each sample. Define a compact kernel function \f$k: \mathbb{R}^{+} \mapsto \mathbb{R}^{+}\f$ such that \f$k(\theta)\f$ is decreasing and \f$k(\theta) = 0\f$ if \f$\theta > 1\f$ (i.e., the support of the kernel is \f$[0,1]\f$). Let
\f{equation*}{
    k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right).
\f}
Define the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is \f$K_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)})\f$.

<B>Configuration Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"NearestNeighbors"   | <tt>YAML::Node</tt> | - | Options for the spi::Tools::NearestNeighbors object.   |
"NumNearestNeighbors"   | <tt>std::size_t</tt> | <tt>10</tt> | The number of nearest neighbors used to compute the bandwidth. |
"ManifoldDimension"   | <tt>double</tt> | <tt>2.0</tt> | The manifold dimension \f$m\f$. |
"TruncationTolerance"   | <tt>double</tt> | If the kernel is a compact kernel (spi::Tools::CompactKernel), the default is \f$\alpha = 1\f$. Otherwise the default is \f$\alpha = -\log{(5 \times 10^{-2})}\f$ | The parameter \f$\alpha\f$ for if we want to truncate the kernel matrix. If the kernel is compact, then this will <em>always</em> be set to one. |
"TruncateKernelMatrix" | <tt>bool</tt> | <tt>true</tt> | <tt>true</tt>: When computing the kernel matrix \f$K_{\epsilon}^{(ij)} = k_{\epsilon}\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right)\f$ only fill the matrix when \f$\frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} < \alpha\f$; <tt>false</tt>: compute and fill every entry in the the kernel matrix (it will be a dense matrix) |
"NumThreads"   | <tt>std::size_t</tt> | <tt>options["NearestNeighbors.NumThreads"]</tt> | The number of <tt>openMP</tt> threads available to this object. |
"BandwidthParameter"   | <tt>double</tt> | <tt>1.0</tt> | The parmeter \f$\epsilon\f$ used to compute the kernel |
"Optimization"   | <tt>YAML::Node</tt> | see options below | Options for the optimization method used to tune algorithm parameters. |

<B>Optimization Parameters:</B>
Parameter Key | Type | Default Value | Description |
------------- | ------------- | ------------- | ------------- |
"Algorithm"   | <tt>std::string</tt> | <tt>"SBPLX"</tt> | Which <a href="https://nlopt.readthedocs.io/en/latest/">NLopt</a> algorithm should we use? Default: <tt>LBFGS</tt> Options: <tt>COBYLA</tt>, <tt>BOBYQA</tt>, <tt>NEWUOA</tt>, <tt>PRAXIS</tt>, <tt>NM</tt>, <tt>SBPLX</tt>, <tt>MMA</tt>, <tt>SLSQP</tt>, <tt>LBFGS</tt>, <tt>PreTN</tt>, <tt>LMVM</tt>  |
"Ftol.AbsoluteTolerance"   | <tt>double</tt> | <tt>1e-6</tt> | Absolute function tolerance.  |
"Ftol.RelativeTolerance"   | <tt>double</tt> | <tt>1e-6</tt> | Relative function tolerance.  |
"Rtol.AbsoluteTolerance"   | <tt>double</tt> | <tt>1e-6</tt> | Absolute state tolerance.  |
"Rtol.RelativeTolerance"   | <tt>double</tt> | <tt>1e-6</tt> | Relative state tolerance.  |
"MaxEvaluations"   | <tt>std::size_t</tt> | <tt>1000</tt> | The maximum number of cost function evaluations.  |
*/
class SampleRepresentation {
public:

  /// Construct the sample representation by sampling a random variable from \f$\psi\f$
  /**
    @param[in] rv The random variable that we wish to sample
    @param[in] options Setup options
  */
  SampleRepresentation(std::shared_ptr<muq::Modeling::RandomVariable> const& rv, YAML::Node const& options);

  /// Construct the sample representation given samples from the underlying distribution \f$\psi\f$
  /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] options Setup options
  */
  SampleRepresentation(std::shared_ptr<muq::SamplingAlgorithms::SampleCollection> const& samples, YAML::Node const& options);

  /// Construct the sample representation given the samples from the underlying distribution \f$\psi\f$
  /**
    @param[in] samples Samples from the underlying distribution \f$\psi\f$
    @param[in] options Setup options
  */
  SampleRepresentation(std::shared_ptr<const spi::Tools::NearestNeighbors> const& samples, YAML::Node const& options);

  virtual ~SampleRepresentation() = default;

  /// Get the bandwith parameter \f$\epsilon\f$
  /**
  \return The bandwith parameter \f$\epsilon\f$
  */
  double BandwidthParameter() const;

  /// How many samples are in the collection?
  /**
    \return The number of samples in the sample collection (SampleRepresentation::samples)
  */
  std::size_t NumSamples() const;

  /// The state dimension
  /**
  \return The state dimension
  */
  std::size_t StateDim() const;

  /// Get the number of nearest neighbors used to compute the bandwidth parameter
  /**
  \return The number of nearest neighbors used to compute the bandwidth parameter
  */
  std::size_t NumNearestNeighbors() const;

  /// Get the kernel \f$k(\theta)\f$ function (spi::Tools::IsotropicKernel)
  /**
  \return A pointer to the kernel \f$k(\theta)\f$ function
  */
  std::shared_ptr<const spi::Tools::IsotropicKernel> Kernel() const;

  /// Get the \f$i^{th}\f$ point \f$\boldsymbol{x}^{(i)}\f$ from the point cloud.
  /**
  @param[in] i We want to get this point
  \return The \f$i^{th}\f$ point from the point cloud.
  */
  Eigen::Ref<Eigen::VectorXd const> Point(std::size_t const i) const;

  /// Construct the (dense) kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  Let \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$ be the squared bandwidth associated with each sample such that \f$I(i,j)\f$ is the index of the \f$j^{th}\f$ furthest index from the point \f$\boldsymbol{x}^{(i)}\f$. By default, this function defines the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is
  \f{equation*}{
    K_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right).
  \f}
  Recall that \f$k(\theta)\f$ is the kernel function.

  This function does NOT compute the kd trees. It assumes that spi::NumericalSolvers::SampleRepresentation::BuildKDTrees has already been called.

  @param[in] eps The parameter \f$\epsilon>0\f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  @param[in] options This parameter is not used in this default implementation, but is used by implementations in its children (defaults to <tt>nullptr</tt>)
  \return Each entry is the sum of a row in the kernel matrix \f$b_i = \sum_{j=1}^{n} K_{\epsilon}^{(ij)}\f$
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat, const void* options = nullptr);

  /// Construct the (dense) kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  This function computes the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is
  \f{equation*}{
    K_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right).
  \f}
  Recall that \f$k(\theta)\f$ is the kernel function.
  @param[in] eps The parameter \f$\epsilon>0\f$
  @param[in] rvec The vector of parameters \f$r_i\f$.
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  \return Each entry is the sum of a row in the kernel matrix \f$b_i = \sum_{j=1}^{n} K_{\epsilon}^{(ij)}\f$
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, Eigen::Ref<Eigen::MatrixXd> kmat);

  /// Construct the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  Let \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$ be the squared bandwidth associated with each sample such that \f$I(i,j)\f$ is the index of the \f$j^{th}\f$ furthest index from the point \f$\boldsymbol{x}^{(i)}\f$. By default, this function defines the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is
  \f{equation*}{
    K_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right).
  \f}
  Recall that \f$k(\theta)\f$ is the kernel function.

  This function does NOT compute the kd trees. It assumes that spi::NumericalSolvers::SampleRepresentation::BuildKDTrees has already been called.

  @param[in] eps The parameter \f$\epsilon>0\f$
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  @param[in] options This parameter is not used in this default implementation, but is used by implementations in its children (defaults to <tt>nullptr</tt>)
  \return Each entry is the sum of a row in the kernel matrix \f$b_i = \sum_{j=1}^{n} K_{\epsilon}^{(ij)}\f$
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat, const void* options = nullptr);

  /// Construct the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  This function computes the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is
  \f{equation*}{
    K_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right).
  \f}
  Recall that \f$k(\theta)\f$ is the kernel function.
  @param[in] eps The parameter \f$\epsilon>0\f$
  @param[in] rvec The vector of parameters \f$r_i\f$.
  @param[out] kmat The kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  \return Each entry is the sum of a row in the kernel matrix \f$b_i = \sum_{j=1}^{n} K_{\epsilon}^{(ij)}\f$
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, Eigen::SparseMatrix<double>& kmat);

  /// Build the kd trees for the nearest neighbor computation
  void BuildKDTrees() const;

  /// Compute the squared bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  /**
  This function does NOT compute the kd trees. It assumes that spi::NumericalSolvers::SampleRepresentation::BuildKDTrees has already been called.

  \return The squared bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  */
  Eigen::VectorXd SquaredBandwidth() const;

  /// Compute the average of the pair-wise kernel derivative evaluations \f$\bar{k}^{\prime} = \frac{d}{d \epsilon} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] \f$
  /**
  We (recursively) compute the average
  \f{equation*}{
  \bar{k}^{\prime} = \frac{d}{d \epsilon} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] = -n^{-2} \sum_{i,j=0}^{n} \left. \frac{d k}{d \theta} \right|_{\theta = \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j}} \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon^2 r_i r_j}
  \f}
  @param[in] eps The bandwidth parameter \f$\epsilon\f$
  @param[in] rvec The vector \f$\boldsymbol{r}\f$
  \return The kernel derivative average \f$\bar{k}^{\prime}\f$.
  */
  virtual double KernelDerivativeAverage(double const eps, Eigen::VectorXd const& rvec) const;

  /// Compute the average of the pair-wise kernel second derivative evaluations \f$\bar{k}^{\prime} = \frac{d^2}{d \epsilon^2} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] \f$
  /**
  We (recursively) compute the average
  \f{equation*}{
  \bar{k}^{\prime} = \frac{d^2}{d \epsilon^2} \left[ n^{-2} \sum_{i,j=0}^{n} k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] = -n^{-2} \sum_{i,j=0}^{n} \left( -2 \left. \frac{d k}{d \theta} \right|_{\theta = \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j}} \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon^3 r_i r_j} - \left. \frac{d^2 k}{d \theta^2} \right|_{\theta = \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j}} \left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon^2 r_i r_j} \right)^2 \right)
  \f}
  @param[in] eps The bandwidth parameter \f$\epsilon\f$
  @param[in] rvec The vector \f$\boldsymbol{r}\f$
  \return The kernel derivative average \f$\bar{k}^{\prime}\f$.
  */
  virtual double KernelSecondDerivativeAverage(double const eps, Eigen::VectorXd const& rvec) const;

  /// Write the samples to file
  /**
    @param[in] filename The name of the file where we are writing the data
    @param[in] dataset The name of the dataset where the samples are stored inside the file (defaults to <tt>"/"</tt>)
  */
  void WriteToFile(std::string const& filename, std::string const& dataset = "/") const;

  /// Get the manifold dimension
  /**
  \return The manifold dimension (could be estimated in child classes)
  */
  double ManifoldDimension() const;

protected:

  /// Compute the entries of the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$
  /**
  Used to construct a sparse kernel matrix

  This function computes the kernel matrix \f$\boldsymbol{K}_{\epsilon}\f$ such that the \f$(i,j)\f$ component is
  \f{equation*}{
    K_{\epsilon}^{(ij)} = k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right).
  \f}
  Recall that \f$k(\theta)\f$ is the kernel function.
  @param[in] eps The parameter \f$\epsilon>0\f$
  @param[in] rvec The vector of parameters \f$r_i\f$.
  @param[out] entries The kernel matrix entries \f$\boldsymbol{K}_{\epsilon}\f$
  \return Each entry is the sum of a row in the kernel matrix \f$b_i = \sum_{j=1}^{n} K_{\epsilon}^{(ij)}\f$
  */
  Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, std::vector<Eigen::Triplet<double> >& entries);

  /// Store the samples from \f$\psi\f$.
  const std::shared_ptr<const spi::Tools::NearestNeighbors> samples;

  /// The number of nearest neighbors used to compute the bandwidth
  const std::size_t numNearestNeighbors;

  /// The kernel function \f$k(\theta)\F$
  /**
  The kernel is defined by a function \f$k: \mathbb{R}^{+} \mapsto \mathbb{R}^{+}\f$ such that \f$k(\theta) = 0\f$ if \f$\theta>1\f$ and
  \f{equation*}{
  k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right),
  \f}
  where \f$\epsilon\f$ and \f$\boldsymbol{r}\f$ are given parameters.
  */
  std::shared_ptr<spi::Tools::IsotropicKernel> kernel;

  /// The dimension of the manifold \f$m\f$
  mutable double manifoldDim;

  /// The number of <tt>openMP</tt> threads available to this object.
  const std::size_t numThreads;

  /// Do we want to truncate the kernel matrix to enforce sparsity
  /**
  <tt>true</tt>: When computing the kernel matrix \f$K_{\epsilon}^{(ij)} = k_{\epsilon}\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right)\f$ only fill the matrix when \f$\frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} < \alpha\f$; <tt>false</tt>: compute and fill every entry in the the kernel matrix (it will be a dense matrix)
  */
  const bool truncateKernelMatrix;

  /// Options for the parameter tuning optimization
  boost::property_tree::ptree pt;

  /// The bandwidth parameter \f$\epsilon\f$
  double bandwidthPara;

private:

  /// Compute the average of the pair-wise kernel derivative evaluations in a given row
  /**
  We compute the average
  \f{equation*}{
  \bar{k}_i^{\prime} = \frac{d}{d \epsilon} \left[ (2(m-k)-(i\in [m,k]))^{-1} \sum_{j=k}^{m} (2.0 - (j==i)) k\left( \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j} \right) \right] = -(2(m-k)-(i\in [m,k]))^{-1} \sum_{j=k}^{m} (2.0 - (j==i)) \left. \frac{d k}{d \theta} \right|_{\theta = \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon r_i r_j}} \frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x}^{(j)}\|^2}{\epsilon^2 r_i r_j}
  \f}
  (note we can also compute the average second dervative---the formula is obtained by differentiating with respect to \f$\epsilon\f$.)
  @param[in] row The \f$i^{th}\f$ row
  @param[in] coli The column \f$k\f$
  @param[in] colj The column \f$m\f$
  @param[in] eps The bandwidth parameter \f$\epsilon\f$
  @param[in] theta The vector \f$\theta_j = \frac{1}{\epsilon r_i r_j}\f$
  @param[in] first <tt>true</tt>: compute the first derivative; <tt>false</tt>: compute the second derivative
  \return The kernel derivative average \f$\bar{k}^{\prime}_i\f$.
  */
  double RecursiveKernelDerivativeRowAverage(std::size_t const row, std::size_t const coli, std::size_t const colj, double const eps, Eigen::VectorXd const& theta, bool const first) const;

  /// Recursively compute the average of the kernel derivative over a given range of rows
  /**
  Compute the sum
  \f{equation*}{
  \bar{k}^{\prime}_{ij} = \sum_{l=i}^{j} \frac{2(n-l)-1)}{(j-i)(2n -j - i)} \bar{k}_l^{\prime}
  \f}
  (see spi::NumericalSolvers::SampleRepresentation::RecursiveKernelDerivativeRowAverage)
  @param[in] rowi The row \f$i\f$
  @param[in] rowj The row \f$j\f$
  @param[in] eps The bandwidth parameter
  @para[in] rvec The variable bandwidth parameters \f$\boldsymbol{r}\f$
  @param[in] first <tt>true</tt>: compute the first derivative; <tt>false</tt>: compute the second derivative
  \return The average kernel derivative
  */
  double KernelDerivativeAverage(std::size_t const rowi, std::size_t const rowj, double const eps, Eigen::VectorXd const& rvec, bool const first) const;

  /// Compute the average of the kernel derivative over a given range of rows
  /**
  Compute the sum
  \f{equation*}{
  \bar{k}^{\prime}_{ij} = \sum_{l=i}^{j} \frac{2(n-l)-1)}{(j-i)(2n -j - i)} \bar{k}_l^{\prime}
  \f}
  (see spi::NumericalSolvers::SampleRepresentation::RecursiveKernelDerivativeRowAverage)
  @param[in] rowi The row \f$i\f$
  @param[in] rowj The row \f$j\f$
  @param[in] eps The bandwidth parameter
  @para[in] rvec The variable bandwidth parameters \f$\boldsymbol{r}\f$
  @param[in] first <tt>true</tt>: compute the first derivative; <tt>false</tt>: compute the second derivative
  \return The average kernel derivative
  */
  double RecursiveKernelDerivativeAverage(std::size_t const rowi, std::size_t const rowj, double const eps, Eigen::VectorXd const& rvec, bool const first) const;

  /// Initialize the sample representation
  /**
  @param[in] options Setup options
  */
  void Initialize(YAML::Node const& options);

  /// Is the kernel a compact kernel?
  bool compactKernel = false;

  /// At what tolerance should we truncate the kernel matrix?
  /**
  Only used if truncateKernelMatrix is <tt>true</tt>.

  When computing the kernel matrix \f$K_{\epsilon}^{(ij)} = k_{\epsilon}\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right)\f$ only fill the matrix when \f$\frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} < \alpha\f$. This is the parameter \f$\alpha\f$.

  If the kernel is a compact kernel (spi::Tools::CompactKernel), the default is \f$\alpha = 1\f$. Otherwise the default is \f$\alpha = -\log{(5 \times 10^{-2})}\f$.

  This choice corresponds to truncting the exponential kernel when
  \f{equation*}{
  \exp{\left(-\frac{\|\boldsymbol{x}^{(i)}-\boldsymbol{x^{(j)}}\|^2}{\epsilon r_i r_j}\right)} \leq 5 \times 10^{-2}
  \f}
  */
  double truncationTol;

  /// The default values for the spi::NumericalSolvers::SampleRepresentation class.
  struct DefaultParameters {
    /// The truncation tolerance for the kernel matrix
    /**
    If the kernel is a compact kernel (spi::Tools::CompactKernel), return \f$1\f$. Otherwise return \f$-\log{(5 \times 10^{-2})}\f$.

    This choice corresponds to truncting the exponential kernel when
    \f{equation*}{
    \exp{\left(-\frac{\|\boldsymbol{x}-\boldsymbol{x^{\prime}}\|^2}{\theta}\right)} \leq 5 \times 10^{-2}
    \f}
    @param[in] compact <tt>true</tt>: This is a compact kernel; <tt>false</tt>: This is not a compact kernel
    \return The default truncation tolerance parameter
    */
    static double TruncationTolerance(bool const compact);

    /// The default bandwidth parameter \f$\epsilon\f$ is \f$1.0\f$
    inline static const double bandwidthPara = 1.0;

    /// The default number of nearest neighbors for the bandwidth computation is \f$10\f$
    inline static const std::size_t numNearestNeighbors = 10;

    /// The default manifold dimension is \f$2\f$
    inline static const double manifoldDim = 2.0;

    /// Default to truncating the kernel matrix
    inline static const bool truncateKernelMatrix = true;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
