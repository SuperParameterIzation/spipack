---
layout: default
---

# Vary collision operator

Consider the collision equation

$$
\begin{equation}
  \partial_t \psi = \phi (P[\psi] - \psi L),
\end{equation}
$$

where

$$
\begin{equation}
  \begin{array}{ccc}
    P[\psi] = \int_{\mathcal{V}} \int_{\mathcal{W}} \omega^{-2} K \psi_{*}^{\prime} \psi_{*} \, d\boldsymbol{W} d \boldsymbol{V^{\prime}} & \mbox{and} & L = \int_{\mathcal{V}} \psi^{\prime} \int_{\mathcal{W}} K \, d\boldsymbol{W} d \boldsymbol{V^{\prime}}
  \end{array}
\end{equation}
$$

are the partial collision operator and the rescaled collision function, respectively. In these examples, we choose constant $K = 1$. We denote $\psi_{\*}^{\prime} = \psi(\boldsymbol{V_{\*}^{\prime}})$, $\psi_{\*} = \psi(\boldsymbol{V_{*}})$, and $\psi^{\prime} = \psi(\boldsymbol{V^{\prime}})$, where the post-collision velocities are

$$
\begin{equation}
  \begin{array}{ccc}
    \boldsymbol{V}_{*} = \boldsymbol{V} + W \boldsymbol{W} & \mbox{and} & \boldsymbol{V}_{*}^{\prime} = \boldsymbol{V}^{\prime} - W \boldsymbol{W}.
  \end{array}
\end{equation}
$$

Here, $\boldsymbol{W}$ is a unit vector samples from $\beta_{\boldsymbol{W}}$ and $W$ is the prescribed post-collision velocity function. We choose $\beta_{\boldsymbol{W}}$ to be a uniform distribution over the unit hypresphere and

$$
\begin{equation}
  W = \frac{1}{2} W_e + \frac{\mbox{sign}(W_e)}{2} \sqrt{ \max{( 0, W_e^2 - 4(1-\gamma) (\boldsymbol{V} \cdot \boldsymbol{V} + \boldsymbol{V^{\prime}} \cdot \boldsymbol{V^{\prime}}) )} },
\end{equation}
$$

where $W_e = -\boldsymbol{W} \cdot (\boldsymbol{V} - \boldsymbol{V^{\prime}})$.
