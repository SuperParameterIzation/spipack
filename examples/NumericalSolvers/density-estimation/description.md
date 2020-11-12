---
layout: default
---

## Density estimation

Let $\psi$ be a probability density function and let $\{ \boldsymbol{x}^{(i)} \}_{i=1}^{n}$ be samples from $\psi$. Define the sample bandwidth

$$
\begin{equation}
  r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2 \tag{Bandwidth}\label{eq:bandwidth},
\end{equation}
$$

where $I(i,j)$ is the index of $j^{th}$ closest sample to $\boldsymbol{x}^{(i)}$. In other words, $r_i$ is the average distance between $\boldsymbol{x}^{(i)}$ and its $k$ nearest neighbors.

For the purpose of this example, assume that $\psi = \mathcal{N}(\boldsymbol{0}, \boldsymbol{I})$ is a 2-dimensional standard normal distribution.

<figure>
<figcaption>The squared bandwidth $r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2$ at all $n=10^4$ samples with $k=25$ nearest neighbors.</figcaption>
<embed src="figures/SquaredBandwidth.pdf" width="500" height="375"
type="application/pdf">
</figure>

Let $k(\theta) = \exp{\left( - \vert \theta \vert \right)}$ and define the kernel

$$
\begin{equation}
k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}) = k\left( \frac{ \| \boldsymbol{x}^{(i)} - \boldsymbol{x}^{(j)} \|^2 }{ \epsilon r_i r_j } \right).
\end{equation}
$$

Let $\epsilon_l \in [\exp{(h l_{min})}, \exp{(h l_{max})}]$ ($h$ is a stepsize parameter) be candidate bandwidth parameters. Define the parameter

$$
\begin{equation}
\Sigma_l = \sum_{i,j=1}^{N} k_{\epsilon}(\boldsymbol{x}^{(i)}, \boldsymbol{x}^{(j)}).
\end{equation}
$$

and an approximation of how it changes with respect to $\epsilon_l$

$$
\begin{equation}
\Sigma_l^{\prime} = \frac{\log{(\Sigma_{l+1})} - \log{(\Sigma_{l})}}{\log{(\epsilon_{l+1})} - \log{(\epsilon_{l})}}
\end{equation}.
$$

Let $\widetilde{\Sigma}_l^{\prime} = \text{max}_{l}({\Sigma_l^{\prime}})$ (with corresponding index $\tilde{l}$). The estimated manifold dimension is $m = 2 \widetilde{\Sigma}_l$; in this case we know that $m = 2$ so we expect $\widetilde{\Sigma}_l = 1$.

<figure>
<figcaption>The estimated density $\psi^{(i)}$ at all $n=10^4$.</figcaption>
<embed src="figures/DensityEstimation.pdf" width="500" height="375"
type="application/pdf">
</figure>

<figure>
<figcaption>The true density evaluation $\mathcal{N}(\boldsymbol{x}^{(i)}; \boldsymbol{0}, \boldsymbol{I})$ at all $n=10^4$.</figcaption>
<embed src="figures/TrueDensity.pdf" width="500" height="375"
type="application/pdf">
</figure>

<embed src="_density-estimation.cpp">
