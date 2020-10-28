#ifndef GRAPHLAPLACIAN_HPP_
#define GRAPHLAPLACIAN_HPP_

namespace spi {
namespace NumericalSolvers {

/// Estimate a weighted Laplacian operator using a graph Laplacian
/**
  Let \f$\psi\f$ be a probability distribution (with corresponding probability density function). Define the weighted Laplacian operator (given some function \f$H\f$)
  \f{equation*}{
    \nabla_{\psi}^2 H = \psi^{-1} \nabla \cdot (\psi \nabla f).
  \f}
  This class uses \f$n\f$ samples \f$\{x^{(i)}\}_{i=1}^{n}\f$ such that \f$x^{(i)} \sim \psi\f$ to approximate the action of this operator and to solve the Poisson equation (given the right hand side \f$R\f$)
  \f{equation*}{
    \nabla_{\psi}^2 H = R
  \f}
  such that \f$\mathbb{E}_{\psi}[H] = 0\f$.
*/
class GraphLaplacian {
public:
  GraphLaplacian() = default;
  virtual ~GraphLaplacian() = default;
private:
};

} // ParticleMethods
} // namespace spi

#endif
