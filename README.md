\section*{Connection between the 2D Hard-Sphere Molecular Dynamics Simulation and the Boltzmann Equation}

\subsection*{1. General Boltzmann Equation}

The Boltzmann equation describes the time evolution of the single-particle distribution function $f(\mathbf{r},\mathbf{v},t)$:
$$
\frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla_\mathbf{r} f + \mathbf{F} \cdot \nabla_\mathbf{v} f = \left(\frac{\partial f}{\partial t}\right)_{\rm coll},
$$
where the left-hand side accounts for the free streaming of particles and the influence of external forces $\mathbf{F}$, and the right-hand side represents the effect of particle collisions.

\subsection*{2. General Form of the Collision Term}

For elastic collisions, the collision term can be written as an integral over all possible collision partners and scattering angles:
$$
\left(\frac{\partial f(\mathbf{v}_1)}{\partial t}\right)_{\rm coll} =
\int \mathrm{d}\mathbf{v}_2 \int \mathrm{d}\Omega \, 
\sigma(\mathbf{v}_1,\mathbf{v}_2,\Omega) \, |\mathbf{v}_1 - \mathbf{v}_2|
\Big[f(\mathbf{v}_1') f(\mathbf{v}_2') - f(\mathbf{v}_1) f(\mathbf{v}_2)\Big],
$$
where:
\begin{itemize}
    \item $\mathbf{v}_1$ and $\mathbf{v}_2$ are the pre-collision velocities,
    \item $\mathbf{v}_1'$ and $\mathbf{v}_2'$ are the post-collision velocities,
    \item $\sigma(\mathbf{v}_1,\mathbf{v}_2,\Omega)$ is the differential cross-section,
    \item $\Omega$ parametrizes the scattering angle.
\end{itemize}
The term $f(\mathbf{v}_1')f(\mathbf{v}_2')$ corresponds to \emph{gain} (particles scattering into $\mathbf{v}_1$), and $f(\mathbf{v}_1)f(\mathbf{v}_2)$ corresponds to \emph{loss} (particles scattering out of $\mathbf{v}_1$).

\subsection*{3. Numerical Approximation in the 2D Hard-Sphere Simulation}

In the molecular dynamics simulation:
\begin{itemize}
    \item $N$ particles are tracked individually with positions $\mathbf{r}_i = (x_i,y_i)$ and velocities $\mathbf{v}_i = (v_{x,i},v_{y,i})$.
    \item Time is discretized with a step $\Delta t$, and particle positions are updated as:
    $$
        \mathbf{r}_i(t+\Delta t) = \mathbf{r}_i(t) + \mathbf{v}_i(t) \Delta t.
    $$
    \item Particle-wall collisions are handled by reflecting velocities, enforcing boundary conditions.
\end{itemize}

\paragraph{Collision approximation:}
\begin{itemize}
    \item The integral over all velocities in the collision term is replaced by a discrete sum over \emph{all particle pairs} $(i,j)$.
    \item A collision occurs if the distance between particles is less than twice the radius:
    $$
        |\mathbf{r}_i - \mathbf{r}_j| < 2r.
    $$
    \item The collision normal is defined as
    $$
        \mathbf{n} = \frac{\mathbf{r}_i - \mathbf{r}_j}{|\mathbf{r}_i - \mathbf{r}_j|},
    $$
    and the relative velocity along the normal is
    $$
        v_n = (\mathbf{v}_i - \mathbf{v}_j) \cdot \mathbf{n}.
    $$
    \item Only approaching particles ($v_n < 0$) undergo velocity updates, corresponding to the \emph{gain and loss} terms:
    $$
        \mathbf{v}_i \rightarrow \mathbf{v}_i - v_n \mathbf{n}, \quad
        \mathbf{v}_j \rightarrow \mathbf{v}_j + v_n \mathbf{n}.
    $$
    \item Overlaps are corrected by pushing particles apart:
    $$
        \mathbf{r}_i \rightarrow \mathbf{r}_i + \frac{2r - |\mathbf{r}_i-\mathbf{r}_j|}{2} \mathbf{n}, \quad
        \mathbf{r}_j \rightarrow \mathbf{r}_j - \frac{2r - |\mathbf{r}_i-\mathbf{r}_j|}{2} \mathbf{n}.
    $$
\end{itemize}

\paragraph{Interpretation:}  
This procedure is a **deterministic, discrete-particle approximation** of the Boltzmann collision term. The continuous integrals over velocities and scattering angles are replaced by sums over actual particle pairs and the collision normal $\mathbf{n}$, with the collision probability determined by geometric overlap. Over many time steps and many particles, the velocity distribution evolves similarly to the prediction of the Boltzmann equation.
