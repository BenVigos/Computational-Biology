# Assigment 3

## The Glazier-Graner-Hogeweg (GGH) Effective Energy

Core Hamiltonian equations that describes the physical shape, volume, and interactions of the cells on the grid.

$$
$$\mathcal{H}_{GGH}=\sum_{\vec{i},\vec{j} \in neighbors}J(\tau(\sigma(\vec{i})),\tau(\sigma(\vec{j})))(1-\delta(\sigma(\vec{i}),\sigma(\vec{j}))) + \sum_{\sigma}\lambda_{vol}(\tau)(v(\sigma)-V_{t}(\tau(\sigma)))^{2}$$
\end{align}
$$

Table: Definitions of the parameters used in the ODEs

| Parameter | Definition |
| :--- | :--- |
| $\vec{i}$ and $\vec{j}$ | The locations of the lattice sites (voxels).|
| $\sigma(\vec{i})$ | The specific cell index located at voxel $\vec{i}$.|
| $\tau(\sigma)$ | The associated cell type of cell $\sigma$.|
|$J$ | The contact energy per unit area between two neighboring cells.|
|$\delta(\sigma(\vec{i}),\sigma(\vec{j}))$ | The Kronecker delta function. It ensures that the adhesion energy is only calculated at the cell boundaries (where the cell indices of neighboring pixels are different).|
|$\lambda_{vol}(\tau)$ | The inverse compressibility of cells of type $\tau$.$v(\sigma)$: The total number of lattice sites (the actual volume) currently occupied by cell $\sigma$.|
$V_{t}(\tau(\sigma))$ | The target volume of the cell. Deviation of the actual volume from this target increases the Effective Energy. |
