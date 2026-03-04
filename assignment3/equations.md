# Spatiotemporal Development of Neural Tissue

## Spatial Signalling in the Development and Function of Neural Connections 
doi: [10.1093/cercor/1.3.199](https://doi.org/10.1093/cercor/1.3.199) ([sci-hub.ru](https://sci-hub.ru/10.1093/cercor/1.3.199))

### Synaptic Voltage

$$v_i(t)=k_i\cdot c_{ij}(t)\cdot(s_j(t)-\theta_e)$$

Where $v_i$ is the postsynaptic voltage, $i$ indexes the postsynaptic site, and $j$ indexes the presynaptic axonal terminal. $s_j(t)$ is the probability that terminal $j$ fires where $s_j\in\{0,1\}$, and $\theta_e$ is an excitatory threshold.

### Cell Firing

Probability of firing:

$$s_i=5.5\cdot\sum_j v_j$$

### Diffusion

$$u_i(\mathbf{x}, t)=D\cdot\nabla^2u(\mathbf{x},t)-k\cdot u(\mathbf{x}, t)+p(\mathbf{x},t)$$

### Synaptic Strenght Change

??

## Our Neuron Network

The neuron network consists of neurons (nodes) in a lattice of a grid. They are connected (edges) on the network.

- $S_{x,y}(t)\in\{0,1\}$: The firing state. It equals 1 if a neuron fires, and 0 if it is quiet.
- $U_{x,y}(t)$: The continuous concentration of our plasticity molecule in that specific square.

- $W_{(i,j)\rightarrow(x,y)}(t)$: The strength of the connection from a neuron at $(i,j)$ to our target neuron at $(x,y)$

The voltage is the following:

$$v_{x,y}(t)=\sum_{i,j}W_{(i,j)\rightarrow(x,y)}(t)\cdot(S_{i,j}(t)-\theta_e)$$

$$P(S_{x,y}(t+1)=1)=\frac{v_{x,y}(t)^n}{\theta_a^n+v_{x,y}(t)^n}$$

- $\theta_a$: The half-activation threshold. This is the input current required for the neuron to have exactly a 50% chance of firing.
- $n$: The Hill coefficient. This controls the "steepness" of the probability curve. If n=1, it's a very gradual, leaky response. As n gets larger (e.g., n=4 or 5), the curve becomes a sharp, sigmoidal step, mimicking the sudden opening of voltage-gated sodium channels.

### Spatio-Temporal Plasticity

The state $A_{(i,j)\rightarrow(x,y)}(t)$ is a boolean tracker of the existence of a connection. The weight of a connection is given by:

$$W^\text{actual}_{(i,j)\rightarrow(x,y)}(t+1)=W_{(i,j)\rightarrow(x,y)}(t)+\eta\cdot S_{i,j}(t)\cdot(U_{x,y}(t)-\theta)$$

Where $\eta$ is the learning rate and $\theta$ is the concentration threshold. The weight of a connection is updated as follows:

$$W_{(i,j)\rightarrow(x,y)}(t+1)=W_{(i,j)\rightarrow(x,y)}(t)+\eta\cdot S_{i,j}(t)\cdot(U_{x,y}(t)-\theta)$$



