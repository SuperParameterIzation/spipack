---
layout: default
---

## Heat matrix eigenvalues

Let $$\psi$$ be a probability density function and let $$\{ \boldsymbol{x}^{(i)} \}_{i=1}^{n}$$ be samples from $$\psi$$. Define the sample bandwidth
$$
\begin{equation}
  r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2
\end{equation}
$$
where $$I(i,j)$$ is the index of $$j^{th}$$ closest sample to $$\boldsymbol{x}^{(i)}$$. In other words, $$r_i$$ is the average distance between $$\boldsymbol{x}^{(i)}$$ and its $$k$$ nearest neighbors.

For the purpose of this example, assume that $$\psi = \mathcal{N}(\boldsymbol{0}, \boldsymbol{I})$$ is a 2-dimensional standard normal distribution.

 <figure>
  <figcaption>The squared bandwidth $$r_i^2 = \frac{1}{k} \sum_{j=1}^{k} \| \boldsymbol{x}^{(i)}-\boldsymbol{x}^{(I(i,j))} \|^2$$ at all $$n=10^4$$ samples with $$k=25$$ nearest neighbors.</figcaption>
  <embed src="figures/SquaredBandwidth.pdf" width="500" height="375"
 type="application/pdf">
</figure>

 ![This is the caption\label{mylabel}](figures/SquaredBandwidth.pdf)
See figure \ref{mylabel}.


Also define the kernel
$$
\begin{equation}
  k(\theta) = \begin{cases}
  1 & \\
  0 & \mbox{if } \theta>1
  \end{cases}
\end{equation}
$$

$$
\begin{equation}
  \xi \label{eq:xi}
\end{equation}
$$

$$
\begin{equation}
  \nabla \psi = 0 \tag{abc}\label{eq:one}
\end{equation}
$$
