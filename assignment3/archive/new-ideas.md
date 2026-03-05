---
geometry:
- margin=1in
---
# Simulating Neuralplasticity using Multiscale Modelling

## Idea and Equations

Here is the complete mathematical framework for your multiscale model. This architecture explicitly links the molecular, intracellular, cellular, and tissue scales using continuous and discrete mathematics.

This hybrid approach—combining Reaction-Diffusion PDEs, Cellular Automata (CA), the Cellular Potts Model (CPM), and Ordinary Differential Equations (ODEs)—is widely used in cutting-edge computational biology platforms like CompuCell3D and Morpheus.

### 1. The Extracellular Scale: Chemical Fields (Coupled PDEs)

The continuous environment contains two primary chemical fields: a chemoattractant/nutrient $C(\mathbf{x}, t)$ guiding the neuron, and a Matrix Metalloproteinase (MMP) enzyme $M(\mathbf{x}, t)$ secreted by the neuron to clear a path.

Their dynamics are governed by coupled reaction-diffusion equations:

$$\frac{\partial C}{\partial t} = D_C \nabla^2 C - \gamma_C C - U_{cell}(C, \mathbf{x}_{cell})$$

$$\frac{\partial M}{\partial t} = D_M \nabla^2 M - \gamma_M M + S_{cell}(M, \mathbf{x}_{membrane})$$

* $D_C$ and $D_M$ are the diffusion coefficients.
* $\gamma_C$ and $\gamma_M$ are the natural decay/clearance rates.
* $U_{cell}$ is the local nutrient uptake by the neuron.
* $S_{cell}$ is the localized secretion of MMPs strictly at the neuron's membrane pixels.

### 2. The Extracellular Matrix (ECM) Scale: (Cellular Automata)

The ECM is a static, rigid scaffold modeled on a discrete CA grid. Let $E(\mathbf{x}, t) \in \{0, 1\}$ represent the state of the ECM at a specific grid coordinate, where $1$ is dense, impenetrable matrix and $0$ is degraded/empty space.

The CA grid updates based on the local concentration of the PDE enzyme field $M(\mathbf{x}, t)$ exceeding a degradation threshold $\theta_M$:

$$E(\mathbf{x}, t+\Delta t) = \begin{cases} 0 & \text{if } M(\mathbf{x}, t) > \theta_M \\ E(\mathbf{x}, t) & \text{otherwise} \end{cases}$$

### 3. The Cellular Morphology Scale: Neuron Dynamics (CPM)

The neuron's physical shape and movement are determined by the stochastic minimization of the system's overall energy, or Hamiltonian ($H$). We construct a specialized Hamiltonian that incorporates basic cell biomechanics, environmental restrictions (from the CA grid), and chemical guidance (from the PDE grid).

$$H = \sum_{\langle \mathbf{i}, \mathbf{j} \rangle} J(\tau_{\sigma_{\mathbf{i}}}, \tau_{\sigma_{\mathbf{j}}})(1 - \delta_{\sigma_{\mathbf{i}}, \sigma_{\mathbf{j}}}) + \lambda_V (v_\sigma - V_{target})^2 + \sum_{\mathbf{i}} \kappa E(\mathbf{i}) \delta_{\sigma_{\mathbf{i}}, cell} - \mu \sum_{\text{membrane}} C(\mathbf{i})$$

* **Term 1 (Adhesion):** $J$ is the surface adhesion cost between neighboring pixels $\mathbf{i}$ and $\mathbf{j}$ of different cell types ($\tau$). The Kronecker delta $\delta$ ensures we only calculate at the cell boundary.
* **Term 2 (Volume constraint):** $\lambda_V$ is the cell's elasticity, keeping its actual volume $v_\sigma$ near its biological target $V_{target}$.
* **Term 3 (ECM CA Coupling):** $\kappa$ is a massive energy penalty. If the neuron tries to copy its pixel into a grid square where the ECM is intact ($E(\mathbf{i}) = 1$), the energy spikes, causing the CPM's Monte Carlo algorithm to reject the movement. The cell *must* degrade the ECM first.
* **Term 4 (PDE Chemotaxis Coupling):** $\mu$ is the chemotactic sensitivity. The energy is lowered if the neuron's membrane moves into pixels with a higher concentration of the chemoattractant $C(\mathbf{i})$.

### 4. The Intracellular Metabolic Scale: Energy (ODEs)

Inside the cell, the absorbed nutrients ($U_{cell}$) from the PDE grid are converted into ATP. We model this using an ODE with Michaelis-Menten enzyme kinetics. Let $A(t)$ be the internal ATP concentration:

$$\frac{dA}{dt} = V_{max} \frac{U_{cell}}{K_m + U_{cell}} - k_{consume} A$$

**Upward Coupling to CPM:** Cell motility requires energy. In the CPM, boundaries fluctuate based on a probability function governed by an "effective temperature" ($T$). We directly couple internal ATP to this temperature: $T(t) = \alpha A(t)$. More ATP means more aggressive growth cone exploration.

### 5. The Gene Regulatory Scale: Protein Synthesis (ODEs)

Finally, we need a Gene Regulatory Network (GRN) to control the secretion of MMPs. We use a Hill-equation ODE to simulate the activation of a gene $G_M$ (which codes for the MMP enzyme) in response to the chemoattractant signal $C$ hitting the cell membrane:

$$\frac{dG_M}{dt} = \beta \frac{C_{membrane}^n}{K_D^n + C_{membrane}^n} - \delta_G G_M$$

**Upward Coupling to PDEs:** The protein produced by this gene dictates the secretion rate into the extracellular space. Therefore, the PDE secretion term $S_{cell}$ from Scale 1 becomes a direct function of the GRN: $S_{cell} \propto G_M(t)$.

---

### Literature & Mathematical Precedents (With DOIs)

Building a model of this magnitude requires robust software platforms (like CompuCell3D or Morpheus) that are explicitly built to solve PDEs, CPM Hamiltonians, and ODEs simultaneously. Here are core papers that outline the mathematical precedents for this exact multiscale architecture:

1. **General Framework for CPM + PDEs + ODEs:**
Swat, M. H., et al. (2012). *Multi-Scale Modeling of Tissues Using CompuCell3D*. Methods in Cell Biology.
**DOI:** [10.1016/B978-0-12-388403-9.00013-8](https://doi.org/10.1016/B978-0-12-388403-9.00013-8)
*(This is the definitive guide on coupling intracellular ODE networks to the CPM and extracellular PDEs).*
2. **Coupling Intracellular Signaling (ODEs) to CPM Morphodynamics:**
Hort, E., et al. (2022). *A multiscale computational model of YAP signaling in epithelial fingering behavior*. Biophysical Journal.
**DOI:** [10.1016/j.bpj.2022.04.010](https://doi.org/10.1016/j.bpj.2022.04.010)
*(Demonstrates how to embed an ODE signaling network directly into a CPM using the Morpheus platform).*
3. **CPM + ECM Degradation + MMP Enzyme Secretion:**
Kusuma, C., et al. (2022). *Multiscale Modeling of MT1-MMP-Mediated Cell Migration: Destabilization of Cell-Matrix Adhesion*. bioRxiv.
**DOI:** [10.1101/2022.10.12.511909](https://www.google.com/search?q=https://doi.org/10.1101/2022.10.12.511909)
*(Excellent mathematical breakdown of coupling internal ODEs for MMP production to CPM membrane dynamics and ECM modification).*
4. **Tip-Cell Migration and Chemotaxis (CPM + PDEs):**
Palm, M. M., et al. (2016). *Computational Screening of Tip and Stalk Cell Behavior Proposes a Role for Apelin Signaling in Sprout Progression*. PLOS ONE.
**DOI:** [10.1371/journal.pone.0159478](https://www.google.com/search?q=https://doi.org/10.1371/journal.pone.0159478)
*(While focused on angiogenesis, the math for leading-edge cell migration driven by PDE chemoattractants is identical to axonal growth cone modeling).*

---

Would you like to explore how to set up the initial grid parameters for one of these specific modeling platforms (like CompuCell3D), or would you rather dive deeper into the Monte Carlo acceptance algorithm that makes the CPM cell actually move?

## Discretization and Parameters

Setting the parameters and discretization scales is arguably the most notoriously difficult part of building a multiscale model. You have to bridge the gap between chemical reactions that happen in milliseconds and cellular movements that take hours, all while ensuring the math doesn't "blow up" due to numerical instability.

Here is a realistic, biologically grounded parameter set for a 2D slice of your neuron model, using standard values found in platforms like CompuCell3D or Morpheus.

---

### 1. Discretization Parameters (Space and Time)

You must decouple your spatial grid from your different temporal clocks. The continuous PDEs (chemistry) need to update hundreds of times for every single update of the CPM (biology).

* **Spatial Resolution ($\Delta x$):** **1 micrometer per pixel**.
* *Why:* A typical neuron soma is roughly 10 to 20 micrometers in diameter, and a growth cone filopodium can be 1 to 2 micrometers wide. Setting 1 pixel = 1 micrometer allows you to capture the fine morphological details of the neuron without making the grid prohibitively massive.


* **Cellular Time Step ($\Delta t_{MCS}$):** **1 minute per Monte Carlo Step (MCS)**.
* *Why:* In the Cellular Potts Model, time is measured in MCS. One MCS is defined as a sequence of random boundary fluctuation attempts equal to the total number of pixels in the grid. Axons typically grow at a rate of roughly 1 to 10 micrometers per hour. Setting 1 MCS to 1 minute matches this biological crawling speed nicely.


* **Chemical Time Step ($\Delta t_{PDE}$):** **0.01 seconds**.
* *Why:* Chemical diffusion is incredibly fast compared to cell crawling. If you try to update the PDE solver only once per minute, the simulation will violently crash. You must satisfy the Von Neumann stability limit for a Forward Euler finite difference scheme:


$$\Delta t_{PDE} \leq \frac{(\Delta x)^2}{4D}$$


Because of this, the PDE solver must run "sub-steps." For every 1 MCS (cellular minute), the PDE solver will execute thousands of internal loops.

---

### 2. Extracellular Chemical Parameters (PDEs)

These represent the morphogens (guidance cues) and the Matrix Metalloproteinases (MMPs) diffusing through the tissue.

| Parameter | Symbol | Value | Biological Meaning |
| --- | --- | --- | --- |
| **Diffusion Coefficient** | $D_C, D_M$ | **10 to 50** micrometers squared per second | Typical diffusion rate of proteins in interstitial fluid. |
| **Decay Rate** | $\gamma_C, \gamma_M$ | **0.001** per second | The natural half-life and breakdown rate of the proteins in the tissue. |
| **Degradation Threshold** | $\theta_M$ | **0.5** (normalized concentration) | The localized concentration of MMPs required to trigger the CA grid to change from "intact matrix" to "degraded". |

---

### 3. Cellular Morphology Parameters (CPM)

These govern the physical "squishiness" and behavior of the neuron on the grid.

| Parameter | Symbol | Value | Biological Meaning |
| --- | --- | --- | --- |
| **Target Area (2D Volume)** | $V_{target}$ | **300** pixels | The ideal 2D footprint of the cell (representing a 300 square micrometer slice). |
| **Volume Elasticity** | $\lambda_V$ | **2.0** | How strictly the cell resists shrinking or swelling. Higher values make the cell rigid; lower values make it highly deformable. |
| **Cell-Fluid Adhesion** | $J(cell, fluid)$ | **20** | The energy cost of the cell boundary touching the surrounding empty fluid. |
| **Cell-ECM Adhesion** | $J(cell, ECM)$ | **5** | The energy cost of touching the rigid CA matrix. (A lower value means the cell prefers to stick to the ECM rather than float in fluid). |
| **Chemotactic Strength** | $\mu$ | **50** | How aggressively the growth cone boundaries bias their random fluctuations up the PDE chemical gradient. |
| **Fluctuation Amplitude** | $T$ | **10 to 50** | The "effective temperature." Directly linked to intracellular ATP. Dictates how wildly the cell membrane writhes and explores. |

---

### 4. Intracellular Kinetics (ODEs)

These parameters govern the internal engine of the cell, linking the absorbed nutrients to ATP production and gene expression.

| Parameter | Symbol | Value | Biological Meaning |
| --- | --- | --- | --- |
| **Maximum Reaction Velocity** | $V_{max}$ | **0.5** micromoles per second | The maximum rate at which the cell's internal enzymes can convert absorbed nutrients into ATP. |
| **Michaelis Constant** | $K_m$ | **0.1** micromoles | The substrate concentration at which the reaction rate is half of $V_{max}$. |
| **Gene Activation Threshold** | $K_D$ | **0.2** (normalized) | The concentration of the external guidance cue required to flip the internal gene regulatory network "ON" to produce MMPs. |
| **Hill Coefficient** | $n$ | **4** | Determines the steepness of the gene activation. A value of 4 creates a sharp, switch-like "ON/OFF" biological response. |
