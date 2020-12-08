---
layout: default
---

# Vary external acceleration

Consider the homogenous advection equation

$$
\begin{equation}
  \partial_t \psi + \nabla_{\boldsymbol{v}} \cdot (\psi \boldsymbol{A}) = 0
\end{equation}
$$


## Linear forcing

$$
\begin{equation}
  \boldsymbol{A} = A_0 (\boldsymbol{U}_{f}-\boldsymbol{V})
\end{equation}
$$

## Quadratic forcing

$$
\begin{equation}
  \boldsymbol{A} = A_0 \| \boldsymbol{U}_{f} - \boldsymbol{V} \| (\boldsymbol{U}_{f}-\boldsymbol{V})
\end{equation}
$$

## Component-wise quadratic forcing

$$
\begin{equation}
  \boldsymbol{A} = A_0 \vert \boldsymbol{U}_f - \boldsymbol{V} \vert * (\boldsymbol{U}_{f}-\boldsymbol{V}),
\end{equation}
$$

where $*$ is the component-wise multiplication operator.
