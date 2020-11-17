---
layout: default
---

## Kolmogorov operator

Let $\psi$ be a probability density function and let $\{ \boldsymbol{x}^{(i)} \}_{i=1}^{n}$ be samples from $\psi$. Given constant $c$ and function f, define the Kolmogorov operator

$$
\begin{equation}
  \mathcal{L}_{\psi, c} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi}
\end{equation}
$$

Note the special cases:
- $c=0$: The Laplacian operator $\mathcal{L}_{\psi,0} f = \Delta f$
- $c=1$: The weighted Laplacian operator $\mathcal{L}_{\psi,1} f = \Delta f + \nabla f \cdot \frac{\nabla \psi}{\psi} = \psi^{-1} \nabla \cdot (\psi \nabla f) = \Delta_{\psi} f$

# Density estimation phase

<figure>
<figcaption>The estimated density $\psi_i \approx \psi(\boldsymbol{x}^{(i)})$ at each sample.</figcaption>
<embed src="figures/DensityEstimation.pdf" width="500" height="375"
type="application/pdf">
</figure>
