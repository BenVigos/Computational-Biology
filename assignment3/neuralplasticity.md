# Project Proposal: Spatial Neural Tissue Model with Gene-Regulated Plasticity

## 1. Project Overview & Biological Rationale

**Objective:** To develop a multi-scale spatiotemporal model of neural plasticity, demonstrating how environmental chemical gradients (BDNF) interact with intracellular gene regulatory networks to drive the physical reorganization of a neural network ("sprouting" and "pruning").

**Biological Background:** The brain wires itself based on the "use it or lose it" principle. Active neurons secrete Brain-Derived Neurotrophic Factor (BDNF). When neighboring neurons or dendritic spines detect high levels of BDNF, they activate specific genes that produce structural proteins (e.g., Actin), prompting them to grow toward the signal. Conversely, in BDNF-deprived areas, these genes are turned off, leading to synaptic pruning.

**Why it fits the assignment:** This model perfectly integrates a **Reaction-Diffusion continuous environment** (BDNF spreading) with a **discrete Cell-Based model** (neurons on a 2D grid), coupled through an **intracellular ODE network** (gene expression driven by Hill kinetics).

---

## 2. Mathematical Framework

The model is built on three interconnected modules that run simultaneously at every time step.

### Module A: The Environment (Reaction-Diffusion PDE)

We model the continuous spatial diffusion of BDNF across a 2D tissue grid using Fick's Second Law.

$$\frac{\partial [BDNF]}{\partial t} = D \nabla^2 [BDNF] + S_{active} - k_{deg}[BDNF]$$

* **$D \nabla^2 [BDNF]$**: The spatial diffusion of BDNF across the 2D lattice.
* **$S_{active}$**: The Source term. It is $>0$ only in grid coordinates occupied by highly active neurons (simulating signal origin).
* **$k_{deg}[BDNF]$**: The Sink term. Represents the natural first-order decay of the neurotrophic factor in the extracellular matrix.

### Module B: The Intracellular Gene Network (ODE)

Inside *every single cell* on the grid, a genetic toggle switch calculates the concentration of a Structural Protein ($P$) based on the local BDNF concentration it senses from the environment.

$$\frac{dP}{dt} = \alpha \frac{[BDNF]^n}{K_d^n + [BDNF]^n} - \gamma P$$

* **$\alpha \frac{[BDNF]^n}{K_d^n + [BDNF]^n}$**: An activating Hill function. The local $[BDNF]$ acts as the transcription factor. If $[BDNF]$ crosses the threshold $K_d$, the synthesis of structural protein $P$ spikes.
* **$\gamma P$**: The natural degradation rate of the structural protein.
* Add proBDNF or MMP9 enzyme

### Module C: Spatial Morphology (Cell-Based / Ising Model)

To determine if a neuron grows (sprouting) or shrinks (pruning), we use an energy-minimization approach inspired by the Hamiltonian / Ising model. The interaction strength ($J$) is dynamically controlled by the internal protein $P$.

$$\Delta H = - J(P) \cdot (\text{Environmental Favorability})$$

* **Rule of Sprouting**: If $P$ is high (above a set threshold), the cell seeks to minimize energy by occupying adjacent empty grid spaces facing the BDNF gradient .
* **Rule of Pruning**: If $P$ drops near zero, the cell cannot maintain its structure, and peripheral grid nodes belonging to the neuron are freed (reverting to empty space).

---

## 3. Implementation Plan (Milestones)

* **Phase 1: The Chemical Landscape (Grid & PDE)**
* Initialize a 2D NumPy array.
* Implement the finite difference method for the Laplacian to simulate BDNF diffusion.
* *Validation:* Ensure a stationary source creates a stable, circular diffusion gradient over time.


* **Phase 2: Single-Cell ODE Integration**
* Place a static "neuron agent" on the grid.
* Write the Euler method solver for the ODE.
* *Validation:* Plot a time-series graph showing the internal protein $P$ rising only when the BDNF wave reaches the cell's specific coordinate.


* **Phase 3: Cell-Based Spatial Updating**
* Implement the Von Neumann or Moore neighborhood checks.
* Write the logic: `If P > threshold -> occupy neighbor node` and `If P < threshold -> vacate node`.
* *Validation:* Observe the cell agent physically changing shape on the 2D grid, growing toward the BDNF source.



---

## 4. Expected Outcomes & Visuals for the Presentation

The final deliverable will feature a side-by-side visual simulation:

1. **A Heatmap of the PDE:** Showing the dynamic waves of BDNF concentration.
2. **A Cellular Automata Grid:** Showing the neural network physically rewiring itself, breaking old connections in low-BDNF zones, and forming dense new clusters in high-BDNF zones.
3. **Phase Plane / Time Series Plots:** Showing the internal ODE dynamics of a selected single cell as it makes the "decision" to grow or shrink.

## References

*   **Solinas, S. M. G., Edelmann, E., Leßmann, V., & Migliore, M. (2019).** A kinetic model for Brain-Derived Neurotrophic Factor mediated spike timing-dependent LTP. *PLoS Computational Biology*, 15(4): e1006975. doi: 10.1371/journal.pcbi.1006975.
*   **Nishiyama, J., & Yasuda, R. (2015).** Biochemical Computation for Spine Structural Plasticity. *Neuron*, 87(1): 63-75. doi: 10.1016/j.neuron.2015.05.043.
*   **Kirchner, J. H., Euler, L., Fritz, I., Ferreira Castro, A., & Gjorgjieva, J. (2025).** Dendritic growth and synaptic organization from activity-independent cues and local activity-dependent plasticity. *eLife*, 12:RP87527. doi: 10.7554/eLife.87527.
*   **Kirchner, J. H., & Gjorgjieva, J. (2021).** Emergence of local and global synaptic organization on cortical dendrites. *Nature Communications*, 12: 4005. doi: 10.1038/s41467-021-23557-3.
