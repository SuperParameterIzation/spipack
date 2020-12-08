---
layout: default
---

# Vary external acceleration

## Linear forcing

$$
\begin{equation}
  \boldsymbol{a} = a_0 (\boldsymbol{u}_{forcing}-\boldsymbol{v})
\end{equation}
$$

## Quadratic forcing

$$
\begin{equation}
  \boldsymbol{a} = a_0 \| \boldsymbol{u} - \boldsymbol{v} \| (\boldsymbol{u}_{forcing}-\boldsymbol{v})
\end{equation}
$$

## Component-wise quadratic forcing

$$
\begin{equation}
  \boldsymbol{a} = a_0 \vert \boldsymbol{u} - \boldsymbol{v} \vert * (\boldsymbol{u}_{forcing}-\boldsymbol{v}),
\end{equation}
$$

where $*$ is the component-wise multiplication operator.
