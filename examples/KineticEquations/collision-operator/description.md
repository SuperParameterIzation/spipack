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

Here, $\boldsymbol{W}$ is a unit vector samples from $\beta_{\boldsymbol{W}}$, $W$ is the prescribed post-collision velocity function, and $\omega = 1 + \boldsymbol{W} \cdot \nabla_{\boldsymbol{V}} W$. We choose $\beta_{\boldsymbol{W}}$ to be a uniform distribution over the unit hypresphere and

$$
\begin{equation}
  W = \frac{1}{2} W_e + \frac{\mbox{sign}(W_e)}{2} \sqrt{ \max{( 0, W_e^2 - 4(1-\gamma) (\boldsymbol{V} \cdot \boldsymbol{V} + \boldsymbol{V^{\prime}} \cdot \boldsymbol{V^{\prime}}) )} },
\end{equation}
$$

where $W_e = -\boldsymbol{W} \cdot (\boldsymbol{V} - \boldsymbol{V^{\prime}})$.

Define the expected kinetic energy

$$
\begin{equation}
  e_{\psi} = \int_{\mathcal{V}} \boldsymbol{V} \cdot \boldsymbol{V} \psi \, d \boldsymbol{V}.
\end{equation}
$$

## Elastic collision

Choose $\gamma = 1$, which implies energy conservation ($\partial_t e_{\psi} = 0$).

<iframe width="1076" height="704" src="https://www.youtube.com/embed/gC_V_I9EuhE" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<figure>
<figcaption>The expected energy $e_{\psi}$.</figcaption>
<embed src="figures/elastic-collisions/ExpectedEnergy.pdf" width="500" height="375"
type="application/pdf">
</figure>

## Inelastic collision

Choose $\gamma = 0.5$, which implies that we lose kinetic energy with each collision.

<iframe width="1076" height="704" src="https://www.youtube.com/embed/MWibkm-fnnU" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<figure>
<figcaption>The expected energy $e_{\psi}$.</figcaption>
<embed src="figures/inelastic-collisions/ExpectedEnergy.pdf" width="500" height="375"
type="application/pdf">
</figure>

## Hyperelastic collision

Choose $\gamma = 1.5$, which implies that we gain kinetic energy with each collision.

<iframe width="1076" height="704" src="https://www.youtube.com/embed/Sha36PDPED4" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<figure>
<figcaption>The expected energy $e_{\psi}$.</figcaption>
<embed src="figures/hyperelastic-collisions/ExpectedEnergy.pdf" width="500" height="375"
type="application/pdf">
</figure>
