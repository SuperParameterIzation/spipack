---
layout: default
---

# Non-zero source term

Consider the evolution equation

$$
\begin{equation}
  \partial_t \psi = \psi R.
\end{equation}
$$

First, define another function $H$ so that

$$
\begin{equation}
  \psi^{-1} \nabla_{\boldsymbol{V}} \cdot (\psi \nabla H) = R.
\end{equation}
$$

We can solve this equation for $H$ and its gradient $\nabla_{\boldsymbol{V}} H$ (see the [Kolmogorov inverse example](../../NumericalSolvers/Kolmogorov-inverse/description.md)). We, therefore, can rewrite the evolution equation

$$
\begin{equation}
  \partial_t \psi - \nabla_{\boldsymbol{V}} \cdot (\psi \nabla H) = 0
\end{equation}
$$

and solve this equation using the same method as the [homogeneous acceleration example](../vary-external-acceleration/description.md).
