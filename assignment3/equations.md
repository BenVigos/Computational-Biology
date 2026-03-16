# Assignment 3

## Simulating Tumor Growth and Angiogenesis in CompuCell3D

## 1. Cellular Potts Model Effective Energy

Core Hamiltonian equation that describes the physical shape, volume, and interactions of the cells on the grid.


$$\mathcal{H}=\sum_{\langle i,j \rangle} J\big(\tau(\sigma_i),\tau(\sigma_j)\big)\left(1-\delta_{\sigma_i,\sigma_j}\right) + \sum_{\sigma} \lambda_V\big(V(\sigma) V_t(\sigma)\big)^2$$

Where:

| Parameter | Definition |
| :--- | :--- |
| $i,j$ | Neighboring lattice sites |
| $\sigma_i$ | Cell identity occupying lattice site $i$ |
| $\tau(\sigma)$ | Type of cell $\sigma$ |
| $J$ | Contact energy between neighboring cell types |
| $\delta$ | Kronecker delta |
| $V(\sigma)$ | Actual cell volume |
| $V_t(\sigma)$ | Target cell volume |
| $\lambda_V$ | Volume constraint strength |

## 3. Diffusion-Reaction Fields

### 3.1 Oxygen field $O$

The oxygen field follows a diffusion-decay equation with type-specific uptake.

$$
\frac{\partial O}{\partial t} = D_O\nabla^2 O - \lambda_O O - \sum_{\tau} U_{\tau}(O)\,\delta_{\tau}
$$

| Parameter | Definition |
| :--- | :--- |
| $D_O$ | Diffusion coefficient for oxygen |
| $\lambda_O$ | Global decay rate of oxygen |
| $U_{\tau}(O)$ | Type-specific oxygen uptake by cells of type $\tau$ |
| $\delta_{\tau}$ | Indicator function equal to 1 at voxels occupied by cell type $\tau$, 0 elsewhere |

Vessel cells act as oxygen reservoirs via fixed-concentration Dirichlet conditions rather than a source term in the PDE:

$$
O = O_{\text{vessel}} \quad \text{at Vascular voxels}
$$

$$
O = O_{\text{neo}} \quad \text{at ActiveNeovascular and InactiveNeovascular voxels}
$$

### 3.2 VEGF1 field $V_1$

VEGF1 is the field secreted by vascular-like cells

$$
\frac{\partial V_1}{\partial t}=D_{V1}\nabla^2V_1-\lambda_{V1}V_1+s_{V1}\,\delta_{\text{vascular-like}}
$$


| Parameter | Definition |
| :--- | :--- |
| $D_{V1}$ | Diffusion coefficient for VEGF1 |
| $\lambda_{V1}$ | Degradation (decay) rate of VEGF1 |
| $s_{V1}$ | Secretion rate of VEGF1 per unit area of source cells |
| $\delta_{\text{vascular-like}}$ | Indicator function equal to 1 at voxels of Vascular, ActiveNeovascular, and InactiveNeovascular cells |

### 3.3 VEGF2 field $V_2$

VEGF2 is the hypoxia-associated signal secreted by hypoxic tumor cells.

$$
\frac{\partial V_2}{\partial t}=D_{V2}\nabla^2V_2-\lambda_{V2}V_2+s_{V2}\,\delta_{\text{Hypoxic}}
$$

| Parameter | Definition |
| :--- | :--- |
| $D_{V2}$ | Diffusion coefficient for VEGF2 |
| $\lambda_{V2}$ | Degradation (decay) rate of VEGF2 |
| $s_{V2}$ | Secretion rate of VEGF2 per unit area of Hypoxic cells |
| $\delta_{\text{Hypoxic}}$ | Indicator function equal to 1 at voxels of Hypoxic cells, 0 elsewhere |

## 4. Tumor Phenotype Switching

Tumor phenotype changes are threshold rules driven by local oxygen at the cell center of mass (COM).

### 4.1 Normal to Hypoxic

$$
O(\mathbf{x}_{\text{COM}}) < O_{\text{nutrient}} \Rightarrow \text{Normal} \to \text{Hypoxic}
$$

### 4.2 Normal or Hypoxic to Necrotic

$$
O(\mathbf{x}_{\text{COM}}) < O_{\text{necrotic}} \Rightarrow \text{Normal/Hypoxic} \to \text{Necrotic}
$$

### 4.3 Hypoxic to Normal recovery

$$
O(\mathbf{x}_{\text{COM}}) > O_{\text{nutrient}} \Rightarrow \text{Hypoxic} \to \text{Normal}
$$

## 5. Tumor Growth Laws

Tumor growth is implemented as a change in target volume, using Michaelis-Menten-like dependence on oxygen.

### 5.1 Normal tumor growth

$$
\frac{dV_t^{N}}{dt} = G_N\frac{O}{K_N + O}
$$

| Parameter | Definition |
| :--- | :--- |
| $G_N$ | Maximum growth rate of a Normal tumor cell |
| $K_N$ | Michaelis-Menten half-saturation constant |
| $O$ | Oxygen concentration |

### 5.2 Hypoxic tumor growth

$$
\frac{dV_t^{H}}{dt} = G_H\frac{O}{K_H + O}
$$

| Parameter | Definition |
| :--- | :--- |
| $G_H$ | Maximum growth rate of a Hypoxic tumor cell |
| $K_H$ | Half-saturation constant for hypoxic growth |
| $O$ | Oxygen concentration |

### 5.3 Necrotic behavior

Necrotic cells do not grow. Their target volume is forced to zero:

$$
V_t^{\text{Necrotic}} = 0
$$


## 6. Vascular Activation and Deactivation

Neovascular sprouts switch phenotype according to the effective VEGF2 level.

### 6.1 Activation

$$
V_{\text{eff}} > \theta_{\text{on}} \Rightarrow \text{InactiveNeovascular} \to \text{ActiveNeovascular}
$$

### 6.2 Deactivation

$$
V_{\text{eff}} < \theta_{\text{off}} \Rightarrow \text{ActiveNeovascular} \to \text{InactiveNeovascular}
$$

## 7. Vascular Growth Law with Contact Inhibition

Vascular target-volume growth is Michaelis-Menten-like.

$$
\frac{dV_t^{\text{neo}}}{dt} = G_V\frac{V_{\text{eff}}}{K_V + V_{\text{eff}}}
$$

| Parameter | Definition |
| :--- | :--- |
| $G_V$ | Maximum growth rate of a neovascular cell |
| $K_V$ | Half-saturation constant |
| $V_{\text{eff}}$ | Effective VEGF2 signal |

But growth only occurs when this condition holds (contact-inhibited sprout growth):

$$
A_{\text{neo-neighbor}} < A_{\text{limit}}
$$

| Parameter | Definition |
| :--- | :--- |
| $A_{\text{neo-neighbor}}$ | Total shared boundary area between the cell and its neovascular (Active or Inactive) neighbors |
| $A_{\text{limit}}$ | Contact-inhibition threshold |


## 8. Mitosis Rules

Cells divide once their actual volume exceeds a type-specific doubling threshold.

### 8.1 Tumor-cell mitosis

$$
V > V_{\text{double,tumor}}
$$

for Normal and Hypoxic cells.

### 8.2 Neovascular mitosis

$$
V > V_{\text{double,vascular}}
$$

for ActiveNeovascular and InactiveNeovascular cells.

## 9. HIF-1α Intracellular Regulatory Network

Local oxygen concentration controls HIF-1α accumulation, and HIF-1α in turn drives VEGF transcriptional activity.

### 9.1 Hypoxia signal from oxygen

Local oxygen is converted to a bounded hypoxia input signal via a Michaelis–Menten-type function.

$$
h(O) = \frac{K_O}{K_O + O}
$$

| Parameter | Definition |
| :--- | :--- |
| $h(O)$ | Hypoxia input signal |
| $K_O$ | Half-saturation constant for oxygen sensing |


### 9.2 HIF-1α regulatory node

HIF-1α is is stabilised by the hypoxia signal and constitutively degraded, following a first-order regulatory ODE.

$$
\frac{dH}{dt} = k_{\text{stab}}\,h(O) - k_{\text{deg}}\,H
$$

| Parameter | Definition |
| :--- | :--- |
| $H$ | Intracellular HIF-1α concentration $\in [0, H_{\max}]$ |
| $k_{\text{stab}}$ | HIF-1α stabilisation rate (production per unit hypoxia signal) |
| $k_{\text{deg}}$ | Constitutive degradation rate (fraction removed per step) |


### 9.3 VEGF transcriptional output node

HIF-1α drives VEGF transcription via a Hill activation function. Tumor-wide HIF state adds a boost to the signal perceived by endothelial cells.

$$
f(H) = \frac{H^n}{K_H^n + H^n}
$$

VEGF transcriptional activity 
$V_{\text{drive}}$ is linearly scaled between a basal floor and a maximum ceiling by this activation:

$$
V_{\text{drive}} = V_{\text{basal}} 
    + \left(V_{\max} - V_{\text{basal}}\right) f(H)
$$


| Parameter | Definition |
| :--- | :--- |
| $f(H)$ | Hill activation function |
| $K_H$ | Hill half-saturation constant for HIF-1α activation |
| $n$ | Hill coefficient |
| $V_{\text{basal}}$ | Basal VEGF transcriptional activity in the absence of HIF-1α |
| $V_{\max}$ | Maximum VEGF transcriptional activity at full HIF-1α saturation |
| $V_{\text{drive}}$ | Intracellular VEGF transcriptional activity proxy |





