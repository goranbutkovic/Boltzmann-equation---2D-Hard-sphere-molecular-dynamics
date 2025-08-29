# 2D Hard-Sphere Molecular Dynamics

This simulation models $N$ particles moving in a 2D box, undergoing elastic collisions with each other and the walls. It approximates the **Boltzmann collision term** deterministically.

---

## 1. Boltzmann Equation

The time evolution of the single-particle distribution function $f(\mathbf{r},\mathbf{v},t)$:

$$
\frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_\mathbf{r} f + \mathbf{F} \cdot \nabla_\mathbf{v} f = \left(\frac{\partial f}{\partial t}\right)_{\rm coll}
$$


> This is a **deterministic, discrete-particle approximation** of the Boltzmann collision term. Over time, the particle velocity distribution evolves similarly to the predictions of the Boltzmann equation.
