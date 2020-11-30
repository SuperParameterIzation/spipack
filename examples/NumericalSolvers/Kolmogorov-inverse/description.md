---
layout: default
---

# Inverse Kolmogorov operator

For more information on the Kolmogorov operator see [the Kolmogorov eigendecomposition example](../Kolmogorov-eigendecomposition/description.md).

Let $\psi$ be a probability density function and let $\left[ \boldsymbol{x}^{(i)} \right]_{i=1}^{n}$ be samples from $\psi$. Given constant $c$ and function f, define the Kolmogorov operator

$$
\begin{equation}
  \mathcal{L}_{\psi, c} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi}
\end{equation}
$$

Note the special cases:
- $c=0$: The Laplacian operator $\mathcal{L}_{\psi,0} f = \Delta f$
- $c=1$: The weighted Laplacian operator $\mathcal{L}_{\psi,1} f = \Delta f + \nabla f \cdot \frac{\nabla \psi}{\psi} = \psi^{-1} \nabla \cdot (\psi \nabla f)$

The goal of this example is to find a function $h$ such that

$$
\begin{equation}
\mathcal{L}_{\psi, c} h = f
\end{equation}
$$

$\mathbb{E}_{\psi}[h] = 0$---this constraint is necessary because the constant function is on the null space of $\mathcal{L}_{\psi}$.
