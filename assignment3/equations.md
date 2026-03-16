---
geometry:
- margin=1in
---
# Assignment 3: Simulating Tumor Growth and Angiogenesis in CompuCell3D

## 1. Cellular Potts Model Effective Energy

Core Hamiltonian equation that describes the physical shape, volume, and interactions of the cells on the grid *(Slide # 26, from Cell-based modelling 2 lecture)*.


$$\begin{equation}
\mathcal{H}=\sum_{\langle i,j \rangle} J\big(\tau(\sigma_i),\tau(\sigma_j)\big)\left(1-\delta_{\sigma_i,\sigma_j}\right) + \sum_{\sigma} \lambda_V\big(V(\sigma) V_t(\sigma)\big)^2
\end{equation}$$

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

## 2. Chemotaxis in the CPM

Directed endothelial migration is represented as an additional energy contribution in the Cellular Potts update rule (*Slide # 46, from Cell-based modelling 2 lecture)*.

$$\begin{equation}
\Delta H_{\text{chem}} = -\lambda_{\text{chem}}\,\big(S(V_2(\mathbf{x}_t))-S(V_2(\mathbf{x}_s))\big)
\end{equation}$$

with a saturating response function

$$\begin{equation}
S(V)=\frac{V}{K_{\text{sat}}+V}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $\lambda_{\text{chem}}$ | Chemotactic sensitivity|
| $V_2(\mathbf{x})$ | VEGF2 concentration at lattice position $\mathbf{x}$ |
| $S(V)$ | Saturating mapping from concentration to chemotactic drive |
| $K_{\text{sat}}$ | Saturation coefficient controlling response compression at high VEGF2 |

so the effective acceptance energy becomes:

$$\begin{equation}
\Delta H_{\text{eff}} = \Delta H_{\text{CPM}} + \Delta H_{\text{chem}}
\end{equation}$$


## 3. Diffusion-Reaction Fields

### 3.1 Oxygen field $O$

The oxygen field follows a diffusion-decay equation with type-specific uptake *(Slide # 44, from Reaction Diffusion Systems lecture)*.


$$\begin{equation}
\frac{\partial O}{\partial t} = D_O\nabla^2 O - \lambda_O O - \sum_{\tau} U_{\tau}(O)\,\delta_{\tau}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $D_O$ | Diffusion coefficient for oxygen |
| $\lambda_O$ | Global decay rate of oxygen |
| $U_{\tau}(O)$ | Type-specific oxygen uptake by cells of type $\tau$ |
| $\delta_{\tau}$ | Indicator function equal to 1 at voxels occupied by cell type $\tau$, 0 elsewhere |

Vessel cells act as oxygen reservoirs via fixed-concentration Dirichlet conditions rather than a source term in the PDE:

$$\begin{equation}
O = O_{\text{vessel}} \quad \text{at Vascular voxels}
\end{equation}$$

$$\begin{equation}
O = O_{\text{neo}} \quad \text{at ActiveNeovascular and InactiveNeovascular voxels}
\end{equation}$$

### 3.2 VEGF1 field $V_1$

VEGF1 is the field secreted by vascular-like cells *(Slide # 44, from Reaction Diffusion Systems lecture)*.

$$\begin{equation}
\frac{\partial V_1}{\partial t}=D_{V1}\nabla^2V_1-\lambda_{V1}V_1+s_{V1}\,\delta_{\text{vascular-like}}
\end{equation}$$


| Parameter | Definition |
| :--- | :--- |
| $D_{V1}$ | Diffusion coefficient for VEGF1 |
| $\lambda_{V1}$ | Degradation (decay) rate of VEGF1 |
| $s_{V1}$ | Secretion rate of VEGF1 per unit area of source cells |
| $\delta_{\text{vascular-like}}$ | Indicator function equal to 1 at voxels of Vascular, ActiveNeovascular, and InactiveNeovascular cells |

### 3.3 VEGF2 field $V_2$

VEGF2 is the hypoxia-associated signal secreted by hypoxic tumor cells *(Slide # 44, from Reaction Diffusion Systems lecture)*.

$$\begin{equation}
\frac{\partial V_2}{\partial t}=D_{V2}\nabla^2V_2-\lambda_{V2}V_2+s_{V2}\,\delta_{\text{Hypoxic}}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $D_{V2}$ | Diffusion coefficient for VEGF2 |
| $\lambda_{V2}$ | Degradation (decay) rate of VEGF2 |
| $s_{V2}$ | Secretion rate of VEGF2 per unit area of Hypoxic cells |
| $\delta_{\text{Hypoxic}}$ | Indicator function equal to 1 at voxels of Hypoxic cells, 0 elsewhere |

## 4. Tumor Phenotype Switching

Tumor phenotype changes are threshold rules driven by local oxygen at the cell center of mass (COM) *(inspired by Spatiotemporal Models lectures, e.g. Slide # 72)*.

### 4.1 Normal to Hypoxic

$$\begin{equation}
O(\mathbf{x}_{\text{COM}}) < O_{\text{nutrient}} \Rightarrow \text{Normal} \to \text{Hypoxic}
\end{equation}$$

### 4.2 Normal or Hypoxic to Necrotic

$$\begin{equation}
O(\mathbf{x}_{\text{COM}}) < O_{\text{necrotic}} \Rightarrow \text{Normal/Hypoxic} \to \text{Necrotic}
\end{equation}$$

### 4.3 Hypoxic to Normal recovery

$$\begin{equation}
O(\mathbf{x}_{\text{COM}}) > O_{\text{nutrient}} \Rightarrow \text{Hypoxic} \to \text{Normal}
\end{equation}$$

## 5. Tumor Growth Laws

Tumor growth is implemented as a change in target volume, using Michaelis-Menten-like dependence on oxygen *(Slide # 9, Spatiotemporal Models 1 lecture; Slide # 59, Enzyme Kinetics lecture)*

### 5.1 Normal tumor growth

$$\begin{equation}
\frac{dV_t^{N}}{dt} = G_N\frac{O}{K_N + O}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $G_N$ | Maximum growth rate of a Normal tumor cell |
| $K_N$ | Michaelis-Menten half-saturation constant |
| $O$ | Oxygen concentration |

### 5.2 Hypoxic tumor growth

$$\begin{equation}
\frac{dV_t^{H}}{dt} = G_H\frac{O}{K_H + O}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $G_H$ | Maximum growth rate of a Hypoxic tumor cell |
| $K_H$ | Half-saturation constant for hypoxic growth |
| $O$ | Oxygen concentration |

### 5.3 Necrotic behavior

Necrotic cells do not grow. Their target volume is forced to zero *(Slide # 25, from Cell-based modelling 2 lecture)*.

$$\begin{equation}
V_t^{\text{Necrotic}} = 0
\end{equation}$$


## 6. Vascular Activation and Deactivation

Neovascular sprouts switch phenotype according to the effective VEGF2 level *(inspired by Spatiotemporal Models lectures, e.g. Slide # 72)*.

### 6.1 Activation

$$\begin{equation}
V_{\text{eff}} > \theta_{\text{on}} \Rightarrow \text{InactiveNeovascular} \to \text{ActiveNeovascular}
\end{equation}$$

### 6.2 Deactivation

$$\begin{equation}
V_{\text{eff}} < \theta_{\text{off}} \Rightarrow \text{ActiveNeovascular} \to \text{InactiveNeovascular}
\end{equation}$$

## 7. Vascular Growth Law with Contact Inhibition

Vascular target-volume growth is Michaelis-Menten-like *(Slide # 9, Spatiotemporal Models 1 lecture; Slide # 59, Enzyme Kinetics lecture)*.

$$\begin{equation}
\frac{dV_t^{\text{neo}}}{dt} = G_V\frac{V_{\text{eff}}}{K_V + V_{\text{eff}}}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $G_V$ | Maximum growth rate of a neovascular cell |
| $K_V$ | Half-saturation constant |
| $V_{\text{eff}}$ | Effective VEGF2 signal |

But growth only occurs when this condition holds (contact-inhibited sprout growth):

$$\begin{equation}
A_{\text{neo-neighbor}} < A_{\text{limit}}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $A_{\text{neo-neighbor}}$ | Total shared boundary area between the cell and its neovascular (Active or Inactive) neighbors |
| $A_{\text{limit}}$ | Contact-inhibition threshold |


## 8. Mitosis Rules

Cells divide once their actual volume exceeds a type-specific doubling threshold *(Slide # 28, from Cell-based modelling 2 lecture)*.

### 8.1 Tumor-cell mitosis

$$\begin{equation}
V > V_{\text{double,tumor}}
\end{equation}$$

for Normal and Hypoxic cells.

### 8.2 Neovascular mitosis

$$\begin{equation}
V > V_{\text{double,vascular}}
\end{equation}$$

for ActiveNeovascular and InactiveNeovascular cells.

## 9. $\text{HIF-1}\alpha$ Intracellular Regulatory Network

Local oxygen concentration controls $\text{HIF-1}\alpha$ accumulation, and $\text{HIF-1}\alpha$ in turn drives VEGF transcriptional activity *(Slide # 43, from Gene Regulatory Networks lecture)*.

### 9.1 Hypoxia signal from oxygen

Local oxygen is converted to a bounded hypoxia input signal via a Michaelis–Menten-type function. 

$$\begin{equation}
h(O) = \frac{K_O}{K_O + O}
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $h(O)$ | Hypoxia input signal |
| $K_O$ | Half-saturation constant for oxygen sensing |


### 9.2 $\text{HIF-1}\alpha$ regulatory node

$\text{HIF-1}\alpha$ is is stabilised by the hypoxia signal and constitutively degraded, following a first-order regulatory ODE.

$$\begin{equation}
\frac{dH}{dt} = k_{\text{stab}}\,h(O) - k_{\text{deg}}\,H
\end{equation}$$

| Parameter | Definition |
| :--- | :--- |
| $H$ | Intracellular $\text{HIF-1}\alpha$ concentration $\in [0, H_{\max}]$ |
| $k_{\text{stab}}$ | $\text{HIF-1}\alpha$ stabilisation rate (production per unit hypoxia signal) |
| $k_{\text{deg}}$ | Constitutive degradation rate (fraction removed per step) |


### 9.3 VEGF transcriptional output node

$\text{HIF-1}\alpha$ drives VEGF transcription via a Hill activation function. Tumor-wide HIF state adds a boost to the signal perceived by endothelial cells.

$$\begin{equation}
f(H) = \frac{H^n}{K_H^n + H^n}
\end{equation}$$

VEGF transcriptional activity 
$V_{\text{drive}}$ is linearly scaled between a basal floor and a maximum ceiling by this activation:

$$\begin{equation}
V_{\text{drive}} = V_{\text{basal}} 
    + \left(V_{\max} - V_{\text{basal}}\right) f(H)
\end{equation}$$


| Parameter | Definition |
| :--- | :--- |
| $f(H)$ | Hill activation function |
| $K_H$ | Hill half-saturation constant for $\text{HIF-1}\alpha$ activation |
| $n$ | Hill coefficient |
| $V_{\text{basal}}$ | Basal VEGF transcriptional activity in the absence of $\text{HIF-1}\alpha$ |
| $V_{\max}$ | Maximum VEGF transcriptional activity at full $\text{HIF-1}\alpha$ saturation |
| $V_{\text{drive}}$ | Intracellular VEGF transcriptional activity proxy |





