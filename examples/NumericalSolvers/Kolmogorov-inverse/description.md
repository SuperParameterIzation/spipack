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
\begin{array}{ccc}
  \mathcal{L}_{\psi, c} h = f & \mbox{and} & \mathbb{E}_{\psi}[h]=0.
\end{array}
\end{equation}
$$

Setting the expected value of $h$ to zero is necessary because the constant function is in the null space of $\mathcal{L}_{\psi}$.

Given the discrete Kolmogorov operator $\boldsymbol{L}_{\psi}$, the discrete version of the problem is

$$
\begin{equation}
    \begin{array}{ccc}
       \boldsymbol{L}_{\psi} \boldsymbol{h} = \boldsymbol{f} & \mbox{and} & \sum_{i=1}^{n} h^{(i)} = 0.
    \end{array}
\end{equation}
$$

Nationally, $\boldsymbol{h}$ and $\boldsymbol{f}$ are vectors such that the $i^{th}$ components $h^{(i)}$ and $f^{(i)}$ are the function evaluations $h^{(i)} = h(\boldsymbol{x}^{(i)})$ and $f^{(i)} = f(\boldsymbol{x}^{(i)})$ at each sample. Recall from [the Kolmogorov eigendecomposition example](../Kolmogorov-eigendecomposition/description.md) that $\boldsymbol{L} = \boldsymbol{S}^{-1} \boldsymbol{\hat{L}} \boldsymbol{S}$, where $\boldsymbol{S}$ is a diagonal matrix and that $\boldsymbol{\hat{L}} \boldsymbol{Q} = \boldsymbol{Q} \boldsymbol{\Lambda}$ is the eigendecomposition of the symmetric matrix $\boldsymbol{\hat{L}}$.

Define coefficients
