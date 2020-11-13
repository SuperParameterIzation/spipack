#ifndef SAMPLEREPRESENTATION_HPP_
#define SAMPLEREPRESENTATION_HPP_

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
"KernelOptions"   | <tt>YAML::Node</tt> | - | The options for the spi::Tools::CompactKernel object |
"TruncationTolerance"   | <tt>double</tt> | If the kernel is a compact kernel (spi::Tools::CompactKernel), the default is \f$\alpha = 1\f$. Otherwise the default is \f$\alpha = -\log{(5 \times 10^{-2})}\f$ | The parameter \f$\alpha\f$ for if we want to truncate the kernel matrix. If the kernel is compact, then this will <em>always</em> be set to one. |
"TruncateKernelMatrix" | <tt>bool</tt> | <tt>true</tt> | <tt>true</tt>: When computing the kernel matrix \f$K_{\epsilon}^{(ij)} = k_{\epsilon}\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right)\f$ only fill the matrix when \f$\frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} < \alpha\f$; <tt>false</tt>: compute and fill every entry in the the kernel matrix (it will be a dense matrix) |
"NumThreads"   | <tt>std::size_t</tt> | <tt>options["NearestNeighbors.NumThreads"]</tt> | The number of <tt>openMP</tt> threads available to this object. |
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

  /// How many samples are in the collection?
  /**
    \return The number of samples in the sample collection (SampleRepresentation::samples)
  */
  std::size_t NumSamples() const;

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
  \return Each entry is the sum of a row in the kernel matrix \f$b_i = \sum_{j=1}^{n} K_{\epsilon}^{(ij)}\f$
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<Eigen::MatrixXd> kmat) const;

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
  Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, Eigen::Ref<Eigen::MatrixXd> kmat) const;

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
  \return Each entry is the sum of a row in the kernel matrix \f$b_i = \sum_{j=1}^{n} K_{\epsilon}^{(ij)}\f$
  */
  virtual Eigen::VectorXd KernelMatrix(double const eps, Eigen::SparseMatrix<double>& kmat) const;

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
  Eigen::VectorXd KernelMatrix(double const eps, Eigen::Ref<const Eigen::VectorXd> const& rvec, Eigen::SparseMatrix<double>& kmat) const;

  /// Build the kd trees for the nearest neighbor computation
  void BuildKDTrees() const;

  /// Compute the squared bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  /**
  This function does NOT compute the kd trees. It assumes that spi::NumericalSolvers::SampleRepresentation::BuildKDTrees has already been called.

  \return The squared bandwidth \f$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2\f$
  */
  Eigen::VectorXd SquaredBandwidth() const;

  /// Write the samples to file
  /**
    @param[in] filename The name of the file where we are writing the data
    @param[in] dataset The name of the dataset where the samples are stored inside the file (defaults to <tt>"/"</tt>)
  */
  void WriteToFile(std::string const& filename, std::string const& dataset = "/") const;

protected:

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

private:

  /// Initialize the sample representation
  /**
  @param[in] options Setup options
  */
  void Initialize(YAML::Node const& options);

  /// Is the kernel a compact kernel?
  bool compactKernel = false;

  /// Do we want to truncate the kernel matrix to enforce sparsity
  /**
  <tt>true</tt>: When computing the kernel matrix \f$K_{\epsilon}^{(ij)} = k_{\epsilon}\left( \frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} \right)\f$ only fill the matrix when \f$\frac{\| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2}{\epsilon r_i r_j} < \alpha\f$; <tt>false</tt>: compute and fill every entry in the the kernel matrix (it will be a dense matrix)
  */
  const bool truncateKernelMatrix;

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

  /// The number of <tt>openMP</tt> threads available to this object.
  const std::size_t numThreads;

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

    /// The default number of nearest neighbors for the bandwidth computation is \f$10\f$
    inline static const std::size_t numNearestNeighbors = 10;

    /// Default to truncating the kernel matrix
    inline static const bool truncateKernelMatrix = true;
  };

  /// Store the default parameter values
  inline static const DefaultParameters defaults;
};

} // namespace NumericalSolvers
} // namespace spi

#endif
