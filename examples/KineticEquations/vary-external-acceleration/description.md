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

add suppose that we initiall draw samples from a standard Gaussian distribution, $\{\boldsymbol{V}^{(i)}\}_{i=1}^{n}$ such that $\boldsymbol{V}^{(i)} \sim \mathcal{N}(\cdot; \boldsymbol{0}, \boldsymbol{I})$.


## Linear forcing

$$
\begin{equation}
  \boldsymbol{A}(\boldsymbol{V}, \boldsymbol{X}, T) = A_0 (\boldsymbol{U}_{f}(\boldsymbol{X}, T)-\boldsymbol{V})
\end{equation}
$$

<iframe width="1076" height="704" src="https://www.youtube.com/embed/S34nyMxcjWk" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Quadratic forcing

$$
\begin{equation}
  \boldsymbol{A}(\boldsymbol{V}, \boldsymbol{X}, T) = A_0 \| \boldsymbol{U}_{f}(\boldsymbol{X}, T) - \boldsymbol{V} \| (\boldsymbol{U}_{f}(\boldsymbol{X}, T)-\boldsymbol{V})
\end{equation}
$$

<iframe width="1076" height="704" src="https://www.youtube.com/embed/gIPCYvI8rXw" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

## Component-wise quadratic forcing

$$
\begin{equation}
  \boldsymbol{A}(\boldsymbol{V}, \boldsymbol{X}, T) = A_0 \vert \boldsymbol{U}_f(\boldsymbol{X}, T) - \boldsymbol{V} \vert * (\boldsymbol{U}_{f}(\boldsymbol{X}, T)-\boldsymbol{V}),
\end{equation}
$$

where $*$ is the component-wise multiplication operator.

<iframe width="1076" height="704" src="https://www.youtube.com/embed/KGyxZ9tOpwg" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
