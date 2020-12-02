---
layout: default
---

# Inverse Kolmogorov operator

For more information on the Kolmogorov operator see [the Kolmogorov eigendecomposition example](../Kolmogorov-eigendecomposition/description.md).

Let $\psi$ be a probability density function and let $\left[ \boldsymbol{x}^{(i)} \right]_{i=1}^{n}$ be samples from $\psi$. Given constant $c$ and function f, define the Kolmogorov operator

$$
\begin{equation}
  \mathcal{L}_{\psi, c} f = \Delta f + c \nabla f \cdot \frac{\nabla \psi}{\psi}
\end{equation}
$$

Note the special cases:
- $c=0$: The Laplacian operator $\mathcal{L}_{\psi,0} f = \Delta f$
- $c=1$: The weighted Laplacian operator $\mathcal{L}_{\psi,1} f = \Delta f + \nabla f \cdot \frac{\nabla \psi}{\psi} = \psi^{-1} \nabla \cdot (\psi \nabla f)$

The goal of this example is to find a function $h$ such that

$$
\begin{equation}
\begin{array}{ccc}
  \mathcal{L}_{\psi, c} h = f & \mbox{and} & \mathbb{E}_{\psi}[h]=0.
\end{array}
\end{equation}
$$

Setting the expected value of $h$ to zero is necessary because the constant function is in the null space of $\mathcal{L}_{\psi}$.

Given the discrete Kolmogorov operator $\boldsymbol{L}_{\psi, c}$, the discrete version of the problem is

$$
\begin{equation}
    \begin{array}{ccc}
       \boldsymbol{L}_{\psi, c} \boldsymbol{h} = \boldsymbol{f} & \mbox{and} & \sum_{i=1}^{n} h^{(i)} = 0.
    \end{array}
\end{equation}
$$

Nationally, $\boldsymbol{h}$ and $\boldsymbol{f}$ are vectors such that the $i^{th}$ components $h^{(i)}$ and $f^{(i)}$ are the function evaluations $h^{(i)} = h(\boldsymbol{x}^{(i)})$ and $f^{(i)} = f(\boldsymbol{x}^{(i)})$ at each sample. Recall from [the Kolmogorov eigendecomposition example](../Kolmogorov-eigendecomposition/description.md) the similarity transformation and eigendecomposition

$$
\begin{equation}
    \begin{array}{ccc}
        \boldsymbol{L}_{\psi, c} = \boldsymbol{S}^{-1} \boldsymbol{\hat{L}}_{\psi, c} \boldsymbol{S} & \mbox{and} & \boldsymbol{\hat{L}}_{\psi, c} \boldsymbol{\hat{Q}} = \boldsymbol{\hat{Q}} \boldsymbol{\Lambda}
    \end{array}
\end{equation}
$$

where $\boldsymbol{S}$ is a diagonal matrix.

Define coefficients

$$
\begin{equation}
    \begin{array}{ccc}
        \boldsymbol{\widetilde{h}} = \boldsymbol{\hat{Q}}^T \boldsymbol{S} \boldsymbol{h} & \mbox{and} & \boldsymbol{\widetilde{f}} = \boldsymbol{\hat{Q}}^T \boldsymbol{S} \boldsymbol{f},
    \end{array}
\end{equation}
$$

allowing us to express $\boldsymbol{h}$ and $\boldsymbol{f}$ as an expansion using the eigenvectors $\boldsymbol{\hat{Q}}$ of $\boldsymbol{\hat{L}}_{\psi, c}$ as a basis. Substituting these coefficients and the eigendecomposition into the discretize problem derives

$$
\begin{equation}
    \boldsymbol{\Lambda} \boldsymbol{\widetilde{h}} = \boldsymbol{\widetilde{f}}.
\end{equation}
$$

Given evaluations of the right hand side at each sample, we compute the coefficients $\boldsymbol{\widetilde{f}}$ and then solutions to the problem are defined by the coefficients $\boldsymbol{\widetilde{h}} = \boldsymbol{\Lambda}^{-\dagger} \boldsymbol{\widetilde{f}}$, where $\boldsymbol{\Lambda}^{-\dagger}$ denotes the pseudo-inverse of the diagonal matrix $\boldsymbol{\Lambda}$. We recover the solution at each sample by expanding using the eigenvectors $\boldsymbol{h} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}} \boldsymbol{\widetilde{h}}$ and shifting the solution so that $\sum_{i=1}^{n} h^{(i)} = 0$.

Given the coefficients and the eigendecomposition we must also compute the unweighted gradient $\nabla_{\boldsymbol{x}} h$. Let $v$ be an additional scalar function. The product rule for the weighted Laplace operator is

$$
\begin{equation}
    \Delta_{\psi} (v h) = v \Delta_{\psi} h + h \Delta_{\psi} v + 2 \nabla_{\boldsymbol{x}} h \cdot \nabla_{\boldsymbol{x}} v.
\end{equation}
$$

Let $(\lambda_j, q_j)$ be eigenvalue/eigenfunction pairs such that $\Delta_{\psi} q_j = \lambda_j q_j$. The product rule implies that

$$
\begin{eqnarray}
    \nabla_{\boldsymbol{x}} q_j \cdot \nabla_{\boldsymbol{x}} q_k &=& \frac{1}{2} ( \Delta_{\psi} (q_j q_k) - q_j \Delta_{\psi} q_k - q_k \Delta_{\psi} q_j ) \\
    &=& \frac{1}{2} ( \Delta_{\psi} (q_j q_k) - (\lambda_j + \lambda_k) q_j q_k )
\end{eqnarray}
$$

Define the coefficient
\begin{equation}
    C_{ljk} = \langle q_l, q_j q_k \rangle_{\psi} = \int_{\mathbb{R}^{2}} q_l^{\prime} q_j q_k \psi \, d\boldsymbol{x} \approx (\boldsymbol{q}_l^{-1})^T (\boldsymbol{q}_j * \boldsymbol{q}_k),
\end{equation}
where $*$ is the component-wise product operator between vectors, $\boldsymbol{q}_{j}$ is the $j^{th}$ eigenvector of the discrete Laplace operator $\boldsymbol{L}_{\psi, c}$ (i.e., the $j^{th}$ column of $\boldsymbol{Q} = \boldsymbol{S}^{-1} \boldsymbol{\hat{Q}}$), and $(\boldsymbol{Q}_l^{-1})^T$ is the $l^{th}$ row of $\boldsymbol{Q}^{-1} = \boldsymbol{\hat{Q}}^T \boldsymbol{S}$. This implies that

$$
\begin{equation}
    Q_j Q_k = \sum_{i=1}^{\infty} C_{ljk} Q_l.
\end{equation}
$$
