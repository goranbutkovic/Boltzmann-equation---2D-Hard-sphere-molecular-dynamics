# 2D Hard-Sphere Molecular Dynamics

This simulation models $N$ particles moving in a 2D box, undergoing elastic collisions with each other and the walls. It approximates the **Boltzmann collision term** deterministically.

---

## 1. Boltzmann Equation

The time evolution of the single-particle distribution function $f(\mathbf{r},\mathbf{v},t)$:

$$
\frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_\mathbf{r} f + \mathbf{F} \cdot \nabla_\mathbf{v} f = \left(\frac{\partial f}{\partial t}\right)_{\rm coll}
$$

where the collision term is

$(\partial f(v_1)/\partial t)_{\rm coll} = \int dv_2 \int d\Omega \, 
\sigma(v_1,v_2,\Omega) \, |v_1 - v_2| \, [f(v_1') f(v_2') - f(v_1) f(v_2)]$
---

## 2. Numerical Approximation

In this code:

- Track positions $\mathbf{r}_i = (x_i,y_i)$ and velocities $\mathbf{v}_i = (v_{x,i},v_{y,i})$.
- Update positions each timestep: 
$$
\mathbf{r}_i(t+\Delta t) = \mathbf{r}_i(t) + \mathbf{v}_i(t)\Delta t
$$
- Detect particle-particle collisions if
$$
|\mathbf{r}_i - \mathbf{r}_j| < 2r
$$
- Update velocities along the collision normal $\mathbf{n} = (\mathbf{r}_i-\mathbf{r}_j)/|\mathbf{r}_i-\mathbf{r}_j|$:
$$
v_n = (\mathbf{v}_i - \mathbf{v}_j) \cdot \mathbf{n}, \quad
\mathbf{v}_i \rightarrow \mathbf{v}_i - v_n \mathbf{n}, \quad
\mathbf{v}_j \rightarrow \mathbf{v}_j + v_n \mathbf{n}.
$$
- Push particles apart to avoid overlap:
$$
\mathbf{r}_i \rightarrow \mathbf{r}_i + \frac{2r - |\mathbf{r}_i-\mathbf{r}_j|}{2} \mathbf{n}, \quad
\mathbf{r}_j \rightarrow \mathbf{r}_j - \frac{2r - |\mathbf{r}_i-\mathbf{r}_j|}{2} \mathbf{n}.
$$

> This is a **deterministic, discrete-particle approximation** of the Boltzmann collision term. Over time, the particle velocity distribution evolves similarly to the predictions of the Boltzmann equation.
