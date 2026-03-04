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

### Anna's notes
Trying to use the most topics...

Active neurons secrete the plasticity molecule $U$ (e.g., BDNF), which diffuses in the 2D extracellular space according to Fick's second law. $$\frac{\partial U_{x,y}}{\partial t} = D \nabla^2 U_{x,y} - k_{deg} U_{x,y} + \gamma S_{x,y}$$, where:
* $D \nabla^2 U$: Spatial diffusion on the grid.
* $- k_{deg} U$: Natural degradation (Sink).
* $+ \gamma S$: Secretion (Source). It activates only when the neuron at $(x,y)$ fires ($S=1$).

Instead of $U$ updating the synapse directly, $U$ acts as a transcription factor for a target gene inside each neuron. The gene produces a structural protein $P$.$$\frac{dP_i}{dt} = \alpha \frac{U_i^m}{K_{GRN}^m + U_i^m} - \beta P_i$$
* The Hill Function: The gene turns on (produces $P$) only if the local concentration of $U$ crosses the threshold $K_{GRN}$.
* Slower dynamics: It acts as a biological "time delay" and noise filter. A brief spike in $U$ won't instantly change the network; the signal must persist long enough to accumulate protein $P$.

Now we update the connection weight $W_{ji}$ (from pre-synaptic neuron $j$ to post-synaptic neuron $i$). We combine our Hebbian learning idea with the new protein $P$ and a spatial energy penalty.$$\Delta W_{ji} = \eta \cdot \left[ P_i \cdot S_j \cdot e^{-\lambda d_{ji}} \right] - \delta W_{ji}$$
* $P_i \cdot S_j$: The weight grows ONLY IF the receiving neuron has built up structural proteins ($P_i > 0$) AND the sending neuron is actively firing ($S_j = 1$).
* $e^{-\lambda d_{ji}}$ (The Spatial Constraint): This is the biophysical energy constraint. $d_{ji}$ is the Euclidean distance between the two neurons on the grid. The exponential decay means short connections are easy to form, but long connections are heavily penalized (simulating the high metabolic cost of long axons).
* $- \delta W_{ji}$: Natural synaptic pruning. If the gene turns off ($P=0$), the weight slowly decays to zero.

