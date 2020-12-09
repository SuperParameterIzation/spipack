---
layout: default
---

# Vary external acceleration

Consider the homogenous advection equation

$$
\begin{equation}
  \partial_t \psi + \nabla_{\boldsymbol{V}} \cdot (\psi \boldsymbol{A}) = 0
\end{equation}
$$

add suppose that we initially draw samples from a standard Gaussian distribution, $\\{\boldsymbol{V}\_0^{(i)}\\}_{i=1}^{n}$ such that $\boldsymbol{V}\_0^{(i)} \sim \mathcal{N}(\cdot; \boldsymbol{0}, \boldsymbol{I})$. We, therefore, update the samples according to

$$
\begin{equation}
  \partial_t V_t = \boldsymbol{A}.
\end{equation}
$$

Below, we show how these samples evolve for different choices of $\boldsymbol{A}$. We also plot the expected energy

$$
\begin{equation}
  e_{\psi} = \int_{\mathcal{V}} \boldsymbol{V} \cdot \boldsymbol{V} \psi \, d \boldsymbol{V}
\end{equation}
$$

## Linear forcing

$$
\begin{equation}
  \boldsymbol{A} = -\boldsymbol{V}
\end{equation}
$$

<iframe width="1076" height="704" src="https://www.youtube.com/embed/S34nyMxcjWk" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<figure>
<figcaption>The expected energy $e_{\psi}$.</figcaption>
<embed src="figures/linear-external-acceleration/ExpectedEnergy.pdf" width="500" height="375"
type="application/pdf">
</figure>

## Quadratic forcing

$$
\begin{equation}
  \boldsymbol{A} = - \| \boldsymbol{V} \| \boldsymbol{V}
\end{equation}
$$

<iframe width="1076" height="704" src="https://www.youtube.com/embed/gIPCYvI8rXw" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<figure>
<figcaption>The expected energy $e_{\psi}$.</figcaption>
<embed src="figures/quadratic-external-acceleration/ExpectedEnergy.pdf" width="500" height="375"
type="application/pdf">
</figure>

## Component-wise quadratic forcing

$$
\begin{equation}

  \boldsymbol{A} = - \vert \boldsymbol{V} \vert * \boldsymbol{V},
\end{equation}
$$

where $*$ is the component-wise multiplication operator.

<iframe width="1076" height="704" src="https://www.youtube.com/embed/KGyxZ9tOpwg" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<figure>
<figcaption>The expected energy $e_{\psi}$.</figcaption>
<embed src="figures/componentwise-quadratic-external-acceleration/ExpectedEnergy.pdf" width="500" height="375"
type="application/pdf">
</figure>
