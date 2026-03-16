# Assigment 3

## Simulating Tumor Growth using CompuCell3D

### The Glazier-Graner-Hogeweg (GGH) Effective Energy

The total energy of the cell configuration that governs all morphodynamics.

$$\begin{align}
H_{\text{GGH}} = \underbrace{\sum_{\vec{i},\,\vec{j}\;\text{neighbors}} J\left(\tau(\sigma(\vec{i})),\,\tau(\sigma(\vec{j}))\right)\left(1 - \delta(\sigma(\vec{i}),\sigma(\vec{j}))\right)}_{\text{contact energy}} + \underbrace{\sum_{\sigma} \lambda_{\text{vol}}(\tau) \left(v(\sigma) - V_t(\tau(\sigma))\right)^2}_{\text{volume constraint}} + \underbrace{\sum_{\sigma} \lambda_{\text{sur}}(\tau) \left(s(\sigma) - S_t(\tau(\sigma))\right)^2}_{\text{surface constraint}}
\end{align}$$


Table: Definitions of the parameters used in the GGH formula. 

| Parameter | Definition |
| :--- | :--- |
|$\vec{i}$ and $\vec{j}$ | The locations of the lattice sites (voxels).|
|$\sigma(\vec{i})$ | The specific cell index located at voxel $\vec{i}$.|
|$\tau(\sigma)$ | The associated cell type of cell $\sigma$.|
|$J(\tau_1, \tau_2)$ | The contact energy between two neighboring cells.|
|$\delta(\sigma(\vec{i}),\sigma(\vec{j}))$ | The Kronecker delta function.|
|$\lambda_{vol}(\tau), \lambda_{sur}(\tau)$ | The inverse compressibility and surface stiffness.|
$V_{t}(\tau(\sigma))$ | The target volume of the cell. |
| $v(\sigma)$ | Actual volume of the cell. |
| $S_t(\tau)$ | The target surface for the cell type. |
| $s(\sigma)$ | Actual surface of the cell | 


### The Chemotaxis Energy
Chemotaxis is added to the Effective Energy to represent the cells movement towards higher chemical concentration.

*Saturated Chemotaxis Rule for VEGF2*

$$
\begin{align}
\Delta\mathcal{H}_{chemotaxis}^{VEGF2} = -\lambda \left(\frac{C(\vec{i}_{target})}{s+C(\vec{i}_{target})}-\frac{C(\vec{i}_{source})}{s+C(\vec{i}_{source})}\right)
\end{align}
$$

| Parameter | Definition |
| :--- | :--- |
|$\lambda$|Chemotaxis strength.|
|$C(\vec{i})$| The concentration of the VEGF2 field at a given target or source.|
|$s$| Saturation Coefficient.|

*Linear Chemotaxis Rule for VEGF1*

$$\begin{align}
\Delta\mathcal{H}_{chemotaxis}^{VEGF1} = -\lambda (C(\vec{i}_{target}) - C(\vec{i}_{source})
\end{align}$$


Then, the total Effective Energy used in Monte Carlo steps: 

$$
\begin{align}
\Delta H_{total} = \Delta H_{chemotaxis} - \Delta_{HGGH}
\end{align}
$$

## Tumor Cell Growth
Tumor cells grow depending on local oxygen concentration.
### Volume Growth
$$\begin{align}
V_{target}​\(t+1)=V_{target}​(t)+k_v\frac{O}{​K+O}​
\end{align}
$$
### Surface Growth
$$\begin{align}
S_{target}​\(t+1)=S_{target}​(t)+k_s\frac{O}{​K+O}​
\end{align}
$$
| Parameter | Definition |
| :--- | :--- |
|$O$ | Oxygen concentration.|
|$K$ |Growth saturation coefficient.|
|$k_v, k_s$ | Volume and surface growth rate.|
### Necrotic Cell Decay
Cells turn necrotic if $O_{local} < \theta_{O, thr}$. Necrotic cells shrink over time.

$$\begin{align}
V_\{target}​(t+1)=max(V_{target}​(t)−k_{necrosis}​,0)
\end{align}$$
## Vascular Growth
Endothelial cells respond to VEGF2 concentration.
### Volume Growth
$$\begin{align}
V_{target}​\(t+1)=V_{target}​(t)+k_v^{vasc}\frac{\text{VEGF}}{​K_v+\text{VEGF}}​
\end{align}
$$
### Surface Growth
$$\begin{align}
S_{target}​\(t+1)=S_{target}​(t)+k_s^{vasc}\frac{\text{VEGF}}{​K_v+\text{VEGF}}​
\end{align}
$$
| Parameter | Definition |
| :--- | :--- |
|$\text{VEGF}$ | VEGF2 concentration.|
|$K_v$ |Growth saturation coefficient.|
|$k_v^{vasc}, k_s^{vasc}$ | Volume and surface growth rate.|

## Cell Division
For both tumor and vascular cells, cells divide when their volume exceeds the volume doubling threshold.

$$\begin{align}
V > V_{\text{doubling}}
\end{align}$$


### Partial Differential Equations for Chemical Fields

### 1. Oxygen Field ($P$)

$$\begin{align}
\frac{\partial P}{\partial t} = D_O \nabla^2 P - \epsilon_{tumor} P \cdot \delta_{tumor} + S_{blood}
\end{align}
$$

| Parameter | Definition |
| :--- | :--- |
|$D_O $| The diffusion constant for the oxygen field.|
|$- \epsilon_{tumor} P$| The rate at which the tumor consumes oxygen.|
|$\delta_{tumor}$ | A function that equals to 1 at voxels where there is a tumoral cell, 0 elsewhere. |
|$S_{blood}$| Oxygen partial pressure constantly supplied by the blood vessels | 

###  2. VEGF2 Field ($V_2$)
$$
\begin{align} 
\frac{\partial V_2}{\partial t}= D_{V_2}\nabla^{2}V_2 -\epsilon_{V_2}V_2+\delta_{hypoxic}x_{V_2}V_2
\end{align}
$$
| Parameter | Definition |
| :--- | :--- |
|$D_{V_2}$| The diffusion constant for the VEGF2 field.|
|$\epsilon_{V_2}$| Decay rate for VEGF2.|
|$x_{V_2}$| Secretion rate.|
|$\delta_{hypoxic}$| A function that equals 1 for hypoxic cells and equals 0 elsewhere.|

###  3. VEGF1 Field ($V_1$)
$$
\begin{align} 
\frac{\partial V_1}{\partial t}= D_{V_1}\nabla^{2}V_1 -\epsilon_{V_1}V_1+\delta_{vasc}x_{V_1}V_1
\end{align}
$$
| Parameter | Definition |
| :--- | :--- |
|$D_{V_1}$| The diffusion constant for the VEGF1 field.|
|$\epsilon_{V_1}$| Decay rate for VEGF1.|
|$x_{V_1}$| Secretion rate.|
|$\delta_{vasc}$| A function that equals 1 for endothelial cells and equals 0 elsewhere.|

### The Gene Regulatory Network HIF - 1 $\alpha$ Dynamics
Added to model if tumoral cell becomes hypoxic.  If oxygen drops (hypoxia), the the constant production ($\alpha_H$) allows to accumulate $H$ inside the cell.


$$\begin{align}
\frac{dH}{dt} = \alpha_H - \beta_H \cdot P \cdot H
\end{align}$$

| Parameter | Definition |
| :--- | :--- |
|$\alpha_H$ | Constant production rate at which the cell continuously produces the HIF - 1 $\alpha$ protein.|
|$\beta_H$| The degradation rate constant.|
|$P$| The local oxygen concentration.|
|$H$| The current concentration of HIF - 1 $\alpha$.|


Mention:
- solver we use
- parameter settings => change them and explain what is the result we biologically observe based on what we saw in class









