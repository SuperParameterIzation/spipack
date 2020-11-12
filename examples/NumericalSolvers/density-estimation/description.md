---
layout: default
---

## Density estimation

Let $\psi$ be a probability density function and let $\{ \boldsymbol{x}^{(i)} \}_{i=1}^{n}$ be samples from $\psi$. Define the sample bandwidth

$$
\begin{equation}
  r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2 \tag{Bandwidth}\label{eq:bandwidth}
\end{equation}
$$

where $I(i,j)$ is the index of $j^{th}$ closest sample to $\boldsymbol{x}^{(i)}$. In other words, $r_i$ is the average distance between $\boldsymbol{x}^{(i)}$ and its $k$ nearest neighbors.

For the purpose of this example, assume that $\psi = \mathcal{N}(\boldsymbol{0}, \boldsymbol{I})$ is a 2-dimensional standard normal distribution.

<figure>
<figcaption>The squared bandwidth $r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2$ at all $n=10^4$ samples with $k=25$ nearest neighbors.</figcaption>
<embed src="figures/SquaredBandwidth.pdf" width="500" height="375"
type="application/pdf">
</figure>

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
