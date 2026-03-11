# Assigment 3

## Simulating Tumor Growth using CompuCell3D

### The Glazier-Graner-Hogeweg (GGH) Effective Energy

Core Hamiltonian equation that describes the physical shape, volume, and interactions of the cells on the grid.

$$
\begin{align}
\mathcal{H}_{GGH}=\sum_{\vec{i},\vec{j} \in neighbors}J(\tau(\sigma(\vec{i})),\tau(\sigma(\vec{j})))(1-\delta(\sigma(\vec{i}),\sigma(\vec{j}))) + \sum_{\sigma}\lambda_{vol}(\tau)(v(\sigma)-V_{t}(\tau(\sigma)))^{2}
\end{align}
$$

Table: Definitions of the parameters used in the GGH formula. 

| Parameter | Definition |
| :--- | :--- |
|$\vec{i}$ and $\vec{j}$ | The locations of the lattice sites (voxels).|
|$\sigma(\vec{i})$ | The specific cell index located at voxel $\vec{i}$.|
|$\tau(\sigma)$ | The associated cell type of cell $\sigma$.|
|$J$ | The contact energy per unit area between two neighboring cells.|
|$\delta(\sigma(\vec{i}),\sigma(\vec{j}))$ | The Kronecker delta function. It ensures that the adhesion energy is only calculated at the cell boundaries (where the cell indices of neighboring pixels are different).|
|$\lambda_{vol}(\tau)$ | The inverse compressibility of cells of type $\tau$.$v(\sigma)$: The total number of lattice sites (the actual volume) currently occupied by cell $\sigma$.|
$V_{t}(\tau(\sigma))$ | The target volume of the cell. Deviation of the actual volume from this target increases the Effective Energy. |

### The Chemotaxis Energy
Chemotaxis is added to the Effective Energy to represent the net effect of the cells preferentially forming pseudopods (cellular extensions) in the direction of the steeper VEGF-A gradient.

$$
\begin{align}
\Delta\mathcal{H}_{chemotaxis}=(\mu(\sigma_{target})-\mu(\sigma_{source})) \left(\frac{V(\vec{i}_{target})}{sV_{0}+V(\vec{i}_{target})}-\frac{V(\vec{i}_{source})}{sV_{0}+V(\vec{i}_{source})}\right)
\end{align}
$$

| Parameter | Definition |
| :--- | :--- |
|$\mu$|The degree of chemotactic response of the specific cell.|
|$V(\vec{i})$| The concentration of the VEGF-A field at a given target or source voxel.|
|$V_{0}$| The threshold value of VEGF-A required to activate inactive neovascular cells.|
|$s$| A positive scaling constant that scales the VEGF-A concentration field relative to the neovascular activation threshold.|

Then, the total Effective Energy: 

$$
\begin{align}
\Delta H_{total} = \Delta H_{chemotaxis} - \Delta_{HGGH}
\end{align}
$$
### Partial Differential Equations for Chemical Fields

###  1. VEGF-A Field ($V$)
$$
\begin{align} 
\frac{\partial V}{\partial t}=-\epsilon_{V}V+\delta_{hypoxic}x_{V}V+D_{V}\nabla^{2}V
\end{align}
$$
| Parameter | Definition |
| :--- | :--- |
|$\epsilon_{V}$| The rate at which the VEGF-A field decays.|
|$x_{V}$| The constant normalized rate at which VEGF-A is secreted.|
|$\delta_{hypoxic}$| A function that equals 1 at voxels belonging to hypoxic cells and equals 0 elsewhere.|
|$D_{V}$| The diffusion constant for the VEGF-A field.|

### 2. Oxygen Field ($P$)

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

### Ordinary Differential Equations for Cell Growth
Michealis-Menten form of how cells grow internally based on the external chemical fields.
### 1. Tumor Cell

$$\begin{align}
\frac{dV_{t}(tumor)}{dt}=\frac{G_{m}pO_{2}(\vec{i})}{O_{0}+pO_{2}(\vec{i})}
\end{align}$$

| Parameter | Definition |
| :--- | :--- |
|$dV_{t}(tumor)/dt$| The rate of increase of the target volume of the normal and hypoxic tumor cells.|
|$G_{m}$| The maximum growth rate of the tumor cell.|
|$pO_{2}(\vec{i})$| The partial pressure of oxygen specifically at the center-of-mass of the cell.|
|$O_{0}$| The Michaelis-Menten constant.|

### 2. Neovascular Cell Proliferation

$$\begin{align}
\frac{dV_{t}(neovacular)}{dt}=\frac{G_{v}V(\vec{i})}{n~V_{0}+V(\vec{i})}
\end{align}$$
| Parameter | Definition |
| :--- | :--- |
|$dV_{t}(neovacular)/dt$| The rate of increase of the target volume for active neovascular cells, accounting for contact-inhibited growth.|
|$G_{v}$ | The maximum growth rate constant.|
|$V(\vec{i})$ | The concentration of VEGF-A at the cell's center-of-mass.|
|$n$ | A scaling constant that describes the proportionality of the activation concentration to the concentration at which the cell's growth rate is half of its maximum.|

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







