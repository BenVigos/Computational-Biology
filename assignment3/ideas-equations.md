## Part 1: Mathematical Model and Discretization

### Continuous Governing Equations

**Extracellular Environment (3D Advection-Diffusion-Reaction)**


$$\frac{\partial N}{\partial t} = D_N \nabla^2 N - \mathbf{v} \cdot \nabla N - k_{chem} N W - \sum_{i=1}^{n} U_{N,i} \delta(\mathbf{x}-\mathbf{x}_i) \tag{1}$$

$$\frac{\partial W}{\partial t} = D_W \nabla^2 W - \mathbf{v} \cdot \nabla W - k_{decay} W + \sum_{i=1}^{n} S_{W,i} \delta(\mathbf{x}-\mathbf{x}_i) \tag{2}$$

**Cell-Environment Interface**


$$U_{N,i} = P_N A_i (N_{ext}(\mathbf{x}_i) - N_{int,i}) \tag{3}$$

**Intracellular Metabolism and GRN**


$$\frac{d N_{int,i}}{dt} = \frac{U_{N,i}}{V_i} - \left( V_{max} \frac{N_{int,i}}{K_M + N_{int,i}} \right) \left( \frac{K_I^n}{K_I^n + R_{B,i}^n} \right) \tag{4}$$

$$\frac{d R_{A,i}}{dt} = \phi \left( V_{max} \frac{N_{int,i}}{K_M + N_{int,i}} \right) \left( \frac{K_I^n}{K_I^n + R_{B,i}^n} \right) - \gamma_A R_{A,i} \tag{5}$$

$$\frac{d R_{B,i}}{dt} = k_{stress} W_{int,i} - \gamma_B R_{B,i} \tag{6}$$

**Morphological Scaling (3D Sphere)**


$$V_i = \frac{1}{6\sqrt{\pi}} A_i^{3/2} \tag{7}$$

$$\frac{dA_i}{dt} = \alpha R_{A,i} - \beta A_i \tag{8}$$

---

### Discretized Numerical Equations

**Spatial Environment (Finite Difference, Upwind Scheme for Advection)**


$$N_{i,j,k}^{m+1} = N_{i,j,k}^m + \Delta t \left[ D_N \left( \frac{\Delta^2 N}{\Delta x^2} + \frac{\Delta^2 N}{\Delta y^2} + \frac{\Delta^2 N}{\Delta z^2} \right)_{i,j,k}^m - v_z \frac{N_{i,j,k}^m - N_{i,j,k-1}^m}{\Delta z} - k_{chem} N_{i,j,k}^m W_{i,j,k}^m - \frac{\sum U_{N,cell}}{\Delta V_{grid}} \right] \tag{9}$$

**Internal Metabolites (Forward Euler)**


$$N_{int,i}^{m+1} = N_{int,i}^m + \Delta t \left[ \frac{U_{N,i}^m}{V_i^m} - \text{Metabolism}(N_{int,i}^m, R_{B,i}^m) \right] \tag{10}$$

$$R_{A,i}^{m+1} = R_{A,i}^m + \Delta t \left[ \phi \cdot \text{Metabolism}(N_{int,i}^m, R_{B,i}^m) - \gamma_A R_{A,i}^m \right] \tag{11}$$

**Growth and Volume Update**


$$A_i^{m+1} = A_i^m + \Delta t (\alpha R_{A,i}^m - \beta A_i^m) \tag{12}$$

$$V_i^{m+1} = \frac{1}{6\sqrt{\pi}} (A_i^{m+1})^{3/2} \tag{13}$$

---

## Part 2: Explanation and Assumptions

### Extracellular Physics (Eq. 1, 2, 9)

The environment is a 3D mineral pore where nutrients ($N$) and waste ($W$) undergo three processes: diffusion ($D$), advection (the upward drift $\mathbf{v}$ from hydrothermal buoyancy), and reaction (chemical degradation).

* **Discretization Logic:** We use a 7-point stencil for the 3D Laplacian to simulate 3D diffusion. The advection term uses an "upwind" scheme, meaning the chemical concentration at a point is influenced by the "wind" blowing from the cell below it. This is crucial for simulating the "wake effect" where downstream cells are poisoned by upstream waste.
* **Assumption:** We assume laminar flow (no turbulence) inside the micro-pore.

### Interface and Scaling (Eq. 3, 7, 13)

Nutrient uptake ($U_N$) is proportional to the surface area ($A$). Because we are in 3D, we assume a spherical geometry.

* **The Square-Cube Constraint:** Equation 7 and 13 link area to volume. As the cell grows, volume increases faster than surface area ($V \propto A^{1.5}$). This creates a physical bottleneck: larger cells have less surface area per unit volume to absorb nutrients, naturally limiting their maximum size.
* **Assumption:** The internal environment is "well-mixed," meaning nutrients diffuse instantly once they cross the membrane.

### Metabolism and GRN (Eq. 4, 5, 6, 10, 11)

This is the "brain" of the cell.

* **Enzyme Kinetics:** Equation 4 uses Michaelis-Menten kinetics, assuming ribozymes ($R_A$) saturate when nutrient levels are very high.
* **Hill Repression:** The second bracket in Equation 5 is a Hill function ($n > 1$). It assumes switch-like behavior: when stress ribozyme $R_B$ hits a threshold, it effectively "turns off" the production of metabolic worker $R_A$.
* **Assumption:** Ribozyme $R_B$ is triggered by internal waste $W_{int}$, which we assume is proportional to the external waste levels $W_{ext}$ at the cell's location.

### Growth and Division (Eq. 8, 12)

The cell grows by converting metabolic output ($R_A$) into lipids/membrane ($A$).

* **Discretization:** We update the area at every time step. If $A_i$ exceeds a threshold $A_{crit}$, the simulation spawns a new agent at the same location with half the area.
* **Assumption:** Maintenance costs ($\beta A$) represent the natural degradation of the lipid membrane over time.

Would you like me to define the **CFL stability criteria** values specifically for this 3D grid to prevent the simulation from crashing?

---

To initialize a 3D spatiotemporal model of this complexity, you must balance the chemical supply from the vent with the metabolic demands of the cells. If the nutrients are too low or the flow is too fast, the population will go extinct (wash-out); if the nutrients are too high, the grid will saturate, and spatial dynamics will become irrelevant.

---

## I. Global Domain and Grid

For a typical mineral pore in an alkaline vent (like Lost City), the scale is microscopic.

| Parameter | Suggested Value | Description |
| --- | --- | --- |
| **Domain Size** | $100 \times 100 \times 500\ \mu\text{m}$ | A tall, narrow 3D rectangular prism. |
| **Grid Resolution** | $\Delta x, \Delta y, \Delta z = 2\ \mu\text{m}$ | Grid cells must be slightly larger than initial cells. |
| **Time Step** | $\Delta t = 0.001 \text{ s}$ | Must satisfy CFL conditions for $D_N$ and $v_z$. |
| **Fluid Velocity** | $v_z = 10-50\ \mu \text{m/s}$ | Upward buoyant drift. |

---

## II. Chemical Initial Conditions

The chimney acts as a "bottom-up" bioreactor.

* **Nutrient Source ($N_{source}$):** $1.0 \text{ mM}$ at the bottom boundary ($z = 0$).
* **Initial Background:** Set $N = 0$ and $W = 0$ throughout the rest of the volume at $t=0$.
* **Boundary Conditions:**
* **Bottom:** Dirichlet ($N = 1.0, W = 0$).
* **Top:** Outflow/Neumann ($\frac{\partial N}{\partial z} = 0, \frac{\partial W}{\partial z} = 0$).
* **Sides:** Periodic or No-Flux (Neumann).



---

## III. Agent (Cell) Initialization

We initialize a "founder" colony. To observe competition and "wake effects," do not scatter them randomly; cluster them.

* **Number of Cells:** $n = 10$ to $20$ agents.
* **Location:** Cluster them in a 3D Gaussian cloud centered at $(L/2, W/2, 20\ \mu\text{m})$. This places them near the source but within the fluid flow.
* **Initial Area ($A_0$):** $100\ \mu\text{m}^2$ (Approx. $5.6\ \mu\text{m}$ diameter sphere).
* **Division Threshold ($A_{crit}$):** $200\ \mu\text{m}^2$.

---

## IV. Suggested Parameter Values (Constants)

These values are derived from prebiotic biophysics (e.g., the Szostak and Mansy labs).

### 1. Transport & Diffusion

* **$D_N, D_W$:** $5.0 \times 10^{-10} \text{ m}^2/\text{s}$ (Small molecules like $H_2$ or nucleotides).
* **$P_N$:** $1.0 \times 10^{-6} \text{ m/s}$ (Permeability of fatty acid membranes).

### 2. Metabolism & Kinetics

* **$V_{max}$:** $1.0 \times 10^{-15} \text{ mol/s}$ per cell.
* **$K_M$:** $0.1 \text{ mM}$ (Nutrient concentration at half-max metabolism).
* **$\gamma_A, \gamma_B$:** $0.01 \text{ s}^{-1}$ (Decay/Degradation of RNA).
* **$K_I$:** $1.0 \times 10^{-18} \text{ mol}$ (Internal $R_B$ threshold for repression).
* **$n$ (Hill coeff):** $4.0$ (High cooperativity for sharp switch-like stress response).

### 3. Growth

* **$\alpha$:** $1.0 \times 10^7\ \mu\text{m}^2/\text{mol}$ (Conversion of processed nutrient to membrane).
* **$\beta$:** $0.001 \text{ s}^{-1}$ (Membrane maintenance/lipid loss).

---

## V. Emergent Behaviors to Look For

1. **The "Toxic Wake":** Cells at the very bottom will flourish. However, their waste ($W$) will drift upward ($+z$). You should see cells at $z = 50\ \mu\text{m}$ struggle or enter a repressed state ($R_B$ high) because they are "eating" the exhaust of the founder cells.
2. **Square-Cube Starvation:** As cells approach $A_{crit}$, watch their $N_{int}$ drop. In 3D, the volume increases so fast that the surface area flux may fail to keep up with $V_{max}$, slowing growth just before division.
3. **Wash-out:** If you increase $v_z$ too much, the nutrients pass the cells too quickly to be absorbed, and the colony will shrink and vanish ($\beta A > \alpha R_A$).