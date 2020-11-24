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

We use the same procedure as in the [density estimation example](../density-estimation/description.md) to compute an estimate of $\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})$ at each sample $\boldsymbol{x}^{(i)}$.

<figure>
<figcaption>The estimated density $\psi^{(i)} \approx \psi(\boldsymbol{x}^{(i)})$ at each sample.</figcaption>
<embed src="figures/DensityEstimation.pdf" width="500" height="375"
type="application/pdf">
</figure>

## Discretize the operator

Parameters:
- The operator constant $c=1$ (defines the operator)
- The exponent $\beta$ (user-prescribed, typically we choose $\beta = -0.5$)
- The bandwidth $\epsilon$ (user-prescribed or tuned using the procedure shown below)
- The manifold dimension $d$ (estimated in the density estimation phase, if not know already)

Let $k(\theta) = \exp{\left( - \vert \theta \vert \right)}$ and define the kernel

$$
\begin{equation}
k_{\epsilon,\beta}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{ \| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2 }{ \epsilon \psi^{ \beta}(\boldsymbol{x}^{(i)}) \psi^{ \beta}(\boldsymbol{x}^{(j)}) } \right).
\end{equation}.
$$

which corresponds to the unnormalized kernel matrix $\boldsymbol{\widetilde{K}}$ such that $\widetilde{K}^{(ij)} = k_{\epsilon,\beta}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)})$.

### Parameter tuning

(Following the same procedure in the [density estimation example](../density-estimation/description.md).)

Let $\epsilon \in [\exp{(l_{min})}, \exp{(l_{max})}]$ be candidate bandwidth parameters. Define the parameter

$$
\begin{equation}
\Sigma_l = \sum_{i,j=1}^{N} \widetilde{K}^{(ij)} = k_{\epsilon,\beta}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}).
\end{equation}
$$

and an approximation of how it changes with respect to $\epsilon_l$

$$
\begin{equation}
\Sigma_l^{\prime} = \frac{\log{(\Sigma_{l+1})} - \log{(\Sigma_{l})}}{\log{(\epsilon_{l+1})} - \log{(\epsilon_{l})}}
\end{equation}.
$$

Let $$\widetilde{\Sigma}_l^{\prime}=\max_{l}{(\Sigma_l^{\prime})}$$ (with corresponding optimal bandwidth parameter $\tilde{\epsilon}$).

<figure>
<figcaption>The value of $\Sigma_l^{\prime}$ for candidate bandwidth parameters $\epsilon$; the optimal bandwidth parameter is $\tilde{\epsilon} \approx 0.28$.</figcaption>
<embed src="figures/LogKernelAvgDerivative.pdf" width="500" height="375"
type="application/pdf">
</figure>

We continue to discretize the operator by defining the normalization

$$
\begin{equation}
q_{\epsilon,\beta} (\boldsymbol{x}^{(i)}) = \psi^{-d \beta}(\boldsymbol{x}^{(i)}) \sum_{j=1}^{n} k_{\epsilon,\beta}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)})
\end{equation}
$$

and the augmented kernel

$$
\begin{equation}
k_{\epsilon,\beta,\alpha}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = \frac{k_{\epsilon,\beta}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)})}{q_{\epsilon, \beta}^{ \alpha}(\boldsymbol{x}^{(j)}) q_{\epsilon, \beta}^{ \alpha}(\boldsymbol{x}^{(i)})}
\end{equation}
$$

where $\alpha = 1 + \frac{1}{2} d \beta + \beta - \frac{1}{2} c$ defines the normalization. Finally, define a second normalization

$$
\begin{equation}
q_{\epsilon,\beta,\alpha} (\boldsymbol{x}^{(i)}) = \sum_{j=1}^{n} k_{\epsilon,\beta, \alpha}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)})
\end{equation}
$$

Denote the symmetric kernel matrix \$\boldsymbol{K}$ such that $K^{(ij)} = k_{\epsilon,\beta,\alpha}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)})$ and the diagonal matrices $\boldsymbol{D}$ and $\boldsymbol{P}$ such that $D^{(ii)} = q_{\epsilon,\beta,\alpha} (\boldsymbol{x}^{(i)})$ and $P^{(ii)} = \psi^{\alpha}(\boldsymbol{x}^{(i)})$. The discrete Kolmogorov operator is

$$
\begin{equation}
\mathcal{L}_{\psi,c} \approx \boldsymbol{L} = \epsilon^{-1} \boldsymbol{P}^{-2} (\boldsymbol{D}^{-1} \boldsymbol{K} - \boldsymbol{I}) = \epsilon^{-1} \boldsymbol{P}^{-2} \boldsymbol{\hat{L}}.
\end{equation}
$$

## Eigen decomposition

Let's take the eigendecomposition of $\boldsymbol{\hat{L}} \boldsymbol{\hat{Q}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda}$.

<figure>
<figcaption>The smallest eigenvalues of $\boldsymbol{\hat{L}}$ and their corresponding eigenfunctions.</figcaption>
<embed src="figures/Eigenfunctions_L.pdf" width="500" height="375"
type="application/pdf">
</figure>

This allows us to define the pseudo-inverse $\boldsymbol{\hat{L}}^{-\dagger} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda}^{-\dagger} \boldsymbol{\hat{Q}}^T$. Therefore, we can approximate the action of the discrete Kolmogorov operator and its pseudo-inverse as $\boldsymbol{L} = \epsilon^{-1} \boldsymbol{P}^{-2} \boldsymbol{\hat{L}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda} \boldsymbol{\hat{Q}}^T$ and $\boldsymbol{L}^{-\dagger} = \epsilon \boldsymbol{\hat{L}}^{-\dagger} \boldsymbol{P}^{2} = \epsilon \boldsymbol{\hat{Q}} \boldsymbol{\Lambda}^{-\dagger} \boldsymbol{\hat{Q}}^T \boldsymbol{P}^{2}$, respectively.

<figure>
<figcaption>The smallest eigenvalues of $\boldsymbol{\hat{L}}$ and their corresponding eigenfunctions.</figcaption>
<embed src="figures/AppliedKolmogorovOperator.pdf" width="500" height="375"
type="application/pdf">
</figure>

## References

- ["Variable bandwidth diffusion kernels" by T. Berry & J. Harlim](https://www.sciencedirect.com/science/article/pii/S1063520315000020)
- ["Data-driven spectral decomposition and forecasting of ergodic dynamical systems" by D. Giannakis](https://www.sciencedirect.com/science/article/pii/S1063520317300982)
