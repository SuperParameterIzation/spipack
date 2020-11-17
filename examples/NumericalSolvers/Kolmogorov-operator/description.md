---
layout: default
---

# Kolmogorov operator

This example relies on work shown in ["Variable bandwidth diffusion kernels" by T. Berry & J. Harlim](https://www.sciencedirect.com/science/article/pii/S1063520315000020).

Let $\psi$ be a probability density function and let $\left[ \boldsymbol{x}^{(i)} \right]_{i=1}^{n}$ be samples from $\psi$. Given constant $c$ and function f, define the Kolmogorov operator

$$
\begin{equation}
  \mathcal{L}_{\psi, c} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi}
\end{equation}
$$

Note the special cases:
- $c=0$: The Laplacian operator $\mathcal{L}_{\psi,0} f = \Delta f$
- $c=1$: The weighted Laplacian operator $\mathcal{L}_{\psi,1} f = \Delta f + \nabla f \cdot \frac{\nabla \psi}{\psi} = \psi^{-1} \nabla \cdot (\psi \nabla f)$


## Density estimation phase

See the [density estimation](../density-estimation/description.md) example for details.

<figure>
<figcaption>The estimated density $\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})$ at each sample.</figcaption>
<embed src="figures/DensityEstimation.pdf" width="500" height="375"
type="application/pdf">
</figure>

Importantly: we have an estimate an $\psi^{(i)}$ of the density at each sample.

## Discretize the operator

Parameters:
- The operator constant $c$
- The exponent $\beta$
- The bandwidth $\epsilon$
- The manifold dimension $d$ (estimated in the density estimation phase, if not know already)

Let $k(\theta) = \exp{\left( - \vert \theta \vert \right)}$ and define the kernel

$$
\begin{equation}
k_{\epsilon,\beta}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{ \| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2 }{ \epsilon \psi^{ \beta}(\boldsymbol{x}^{(i)}) \psi^{ \beta}(\boldsymbol{x}^{(j)}) } \right).
\end{equation}.
$$

### Parameter tuning

Following the same procedure in the [density estimation](../density-estimation/description.md) example.

Let $\epsilon \in [\exp{(l_{min})}, \exp{(l_{max})}]$ be candidate bandwidth parameters. Define the parameter

$$
\begin{equation}
\Sigma_l = \sum_{i,j=1}^{N} k_{\epsilon,\beta}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}).
\end{equation}
$$

and an approximation of how it changes with respect to $\epsilon_l$

$$
\begin{equation}
\Sigma_l^{\prime} = \frac{\log{(\Sigma_{l+1})} - \log{(\Sigma_{l})}}{\log{(\epsilon_{l+1})} - \log{(\epsilon_{l})}}
\end{equation}.
$$

Let $$\widetilde{\Sigma}_l^{\prime}=\max_{l}{(\Sigma_l^{\prime})}$$ (with corresponding index $\tilde{l}$ and optimal bandwidth parameter $\tilde{\epsilon}$).

<figure>
<figcaption>The value of $\Sigma_l^{\prime}$ for candidate bandwidth parameters $\epsilon$; the optimal bandwidth parameter is $\tilde{\epsilon} \approx 0.28$.</figcaption>
<embed src="figures/LogKernelAvgDerivative.pdf" width="500" height="375"
type="application/pdf">
</figure>

## References

- ["Variable bandwidth diffusion kernels" by T. Berry & J. Harlim](https://www.sciencedirect.com/science/article/pii/S1063520315000020)
- ["Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis](https://www.sciencedirect.com/science/article/pii/S1063520317300982)
