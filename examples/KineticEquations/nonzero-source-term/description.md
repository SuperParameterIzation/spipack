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
  \psi^{-1} \nabla_{\boldsymbol{V}} \cdot (\psi \nabla_{\boldsymbol{V}} H) = R.
\end{equation}
$$

We can solve this equation for $H$ and its gradient $\nabla_{\boldsymbol{V}} H$ (see the [Kolmogorov inverse example](../../NumericalSolvers/Kolmogorov-inverse/description.md)). We, therefore, can rewrite the evolution equation

$$
\begin{equation}
  \partial_t \psi - \nabla_{\boldsymbol{V}} \cdot (\psi \nabla_{\boldsymbol{V}} H) = 0
\end{equation}
$$

and solve this equation using the same method as the [homogeneous acceleration example](../vary-external-acceleration/description.md). In this example, define $R = -\boldsymbol{V} \cdot \boldsymbol{Y}$, where $\boldsymbol{Y} = [1, 0]$.

<iframe width="1076" height="704" src="https://www.youtube.com/embed/9vUKS9TgS3E" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<figure>
<figcaption>The expected energy $e_{\psi} = \int_{\mathcal{V}} \boldsymbol{V} \cdot \boldsymbol{V} \psi \, d \boldsymbol{V}$.</figcaption>
<embed src="figures/ExpectedEnergy.pdf" width="500" height="375"
type="application/pdf">
</figure>
