---
geometry:
- margin=1in
---
# Assignment 2: Model development and equations 

## Viterbi algorithm
The Viterbi algorithm has the following steps:

1. *Initialization:* 
$$\delta_1(s) = \Pi(s) \times B(O_1|s)$$
where $\Pi(s)$ is the initial probability of being in state s and $B(O_1|s)$ is the emission probability of observing $O_1$ given state s.

and $$\psi_1(s) = none$$

2. *Recursion:* 
$$\delta_t(s) = \max_{r\in s} [\delta_{t-1}(r) \times A(r,s), B(O_t|s)]$$
where $A(r,s)$ is the transition probability from state r to state s.

and $$\psi_t(s) = argmax [\delta_{t}(r) \times A(r,s)]$$
3. *Termination:* 
$$P_t = \max{[\delta_t(s)]}$$
4. *Backtracking:* Trace back through the states to determine the most likely sequence of states.

where above, $\delta_t(s)$ is the maximum probability of being in state s at time t given the observed sequence up to time t,
$B(O_t|s)$ is the emission probability of observing $O_t$ given state s,
and $\psi_t(s)$ is the state at time t-1 that maximizes the probability of being in state s at time t.


## Question 1

$$
\begin{align}
h^+(P_i, \theta_i, n_i) &= \frac{P_i^n}{\theta_i^n + P_i^n} \\
h^-(P_i, \theta_i, n_i) &= 1 - h^+(P_i, \theta_i, n_i) \\
&= \frac{\theta_i^n}{\theta_i^n + P_i^n}
\end{align}
$$

### ODEs - Mechanism I: Transcriptional Hijack

$$
\begin{align}
\frac{dr_a}{dt} &= m_a \times h^+(P_b, \theta_b, n_b) + \gamma_a \times r_a \\
\frac{dr_b}{dt} &= m_b \times h^-(P_a, \theta_a, n_a) + \gamma_b \times r_b \\
\frac{dP_a}{dt} &= k_{Pa} \times r_a - \delta_{Pa} \times P_a \\  
\frac{dP_b}{dt} &= k_{Pb} \times r_b - \delta_{Pb} \times P_b
\end{align}
$$

When transcriptional hijack, $h^-(P_a, \theta_a, n_a) = 1$

Table: Definitions of the parameters used in the ODEs

| Parameter | Definition |
| :--- | :--- |
| $r_i$ | Concentration of transcribed mRNA A or B (M, mol L$^{-1}$). |
| $\dfrac{dr_i}{dt}$ | Time-derivative / change in concentration of transcribed mRNA (M s$^{-1}$). |
| $m_i$ | Maximum transcription rate coefficient for mRNA (M s$^{-1}$). |
| $h^+(P, \theta, n)$ | Hill activation function (dimensionless). |
| $h^-(P, \theta, n)$ | Hill inhibition function (dimensionless). |
| $\theta_i$ | Expression threshold for protein binding (M, mol L$^{-1}$). |
| $n_i$ | Hill coefficient representing regulatory nonlinearity (dimensionless, unitless). |
| $\gamma_i$ | Degradation rate of mRNA (s$^{-1}$). |
| $P_i$ | Concentration of Protein (M, mol L$^{-1}$). |
| $\dfrac{dP_i}{dt}$ | Time-derivative / change in protein concentration (M s$^{-1}$). |
| $k_{Pi}$ | Effective translation rate (s$^{-1}$) converting mRNA concentration to protein production rate. |
| $\delta_{Pi}$ | Degradation rate of Protein A or Protein B (s$^{-1}$). |
| $t$ | Time (s, seconds). |

### SDEVelo - Mechanism II: Splicing Sabotage

$$
\begin{align}
\alpha_A(t) &= \frac{c_A}{1+\exp b_A(t-a_A)} \\
\beta_A^* &= \beta_A h^+(P_B, \theta_B, n_B) \\
dU_A &= (\alpha_A(t) - \beta_A^* U_A(t))dt + \sigma_{1A}dW_{1A} \\
dS_A &= (\beta_A^* U_A(t) - \gamma_A S_A(t))dt + \sigma_{2A}dW_{2A} \\
dP_A &= (k_{PA} S_A(t) - \delta_{PA} P_A(t))dt \\
\end{align}
$$

$$
\begin{align}
\alpha_B(t) &= \frac{c_B}{1+\exp b_B(t-a_B)} \\
\beta_B^* &= \beta_B h^-(P_A, \theta_A, n_A) \\
dU_B &= (\alpha_B(t) - \beta_B^* U_B(t))dt + \sigma_{1B}dW_{1B} \\
dS_B &= (\beta_B^* U_B(t) - \gamma_B S_B(t))dt + \sigma_{2B}dW_{2B} \\
dP_B &= (k_{PB} S_B(t) - \delta_{PB} P_B(t))dt \\
\end{align}
$$

Table: Definitions of the parameters used in the SDEVelo equations

| Parameter | Definition |
| :--- | :--- |
| $\alpha_i(t)$ | Time-dependent transcription rate of unspliced mRNA (pre-mRNA) (M s$^{-1}$). |
| $c_i$ | Maximum transcription rate coefficient for pre-mRNA (M s$^{-1}$). |
| $b_i$ | Steepness parameter of the transcription rate sigmoid function (s$^{-1}$). |
| $t$ | Time (s). |
| $a_i$ | Time shift or activation delay for transcription (s). |
| $\beta^*_i$ | Regulated splicing rate influenced by protein interactions (s$^{-1}$). |
| $\beta_i$ | Base splicing rate parameter (s$^{-1}$). |
| $h^+(P, \theta, n)$ | Hill activation function (dimensionless). |
| $h^-(P, \theta, n)$ | Hill inhibition function (dimensionless). |
| $P_i$ | Concentration of Protein (M). |
| $\theta_i$ | Expression threshold for protein binding (M). |
| $n_i$ | Hill coefficient (dimensionless). |
| $U_i(t)$ | Concentration of unspliced mRNA (pre-mRNA) at time $t$ (M). |
| $dU_i$ | Change in unspliced mRNA concentration (M s$^{-1}$). |
| $S_i(t)$ | Concentration of spliced mRNA at time $t$ (M). |
| $dS_i$ | Change in spliced mRNA concentration (M s$^{-1}$). |
| $\gamma_i$ | Degradation rate of spliced mRNA (s$^{-1}$). |
| $\sigma_{1i}$ | Noise intensity parameter for pre-mRNA transcription (M s$^{-1/2}$). |
| $\sigma_{2i}$ | Noise intensity parameter for the splicing process (M s$^{-1/2}$). |
| $dW_{1i}, dW_{2i}$ | Differentials of the Wiener process (units of s$^{1/2}$); stochastic increments scale as $\sqrt{\Delta t}$. |
| $k_{Pi}$ | Translation rate (s$^{-1}$). |
| $\delta_{Pi}$ | Protein degradation rate (s$^{-1}$). |
| $dt$ | Time differential (s). |


## Downstream metabolic effects
$$
\begin{align}
\frac{dR}{dt} &= \alpha R - \beta RE \\
\frac{dE}{dt} &= -\gamma E + \delta RE
\end{align}
$$

Where:

* $\alpha R$: represents the natural growth or replenishment of the resource in the absence of any interactions with the enzyme E. 
* $-\beta R E$: represents the rate at which the enzyme E utilizes the resource R. 
* $-\gamma E$: represents the rate at which the enzyme E is lost or deactivated over time, independent of its interaction with the resource R. 
* $\delta R E$: represents the rate at which the enzyme E is generated or activated in response to the presence of the resource R. 

with $\alpha=2\;[\mathrm{T^{-1}}]$, $\beta=1.1\;[\mathrm{M^{-1} L^{3} T^{-1}}]$, $\gamma=1\;\mathrm{[T^{-1}]}$, $\delta=0.9\;\mathrm{[M^{-1} L^{3} T^{-1}]}$, $R(0)=1\;\mathrm{[M^{1} L^{-3}]}$ and $E(0)=0.5\;\mathrm{[M^{1} L^{-3}]}$

The fixed points of this system can be found by setting the derivatives to zero:
$$
\begin{align}
0 &= \alpha R - \beta R E \\
0 &= -\gamma E + \delta E R
\end{align}
$$

and solving for R and E. This yields the fixed points:
1. $(R^*, E^*) = (0, 0)$
2. $(R^*, E^*) = \left(\frac{\gamma}{\delta}, \frac{\alpha}{\beta}\right)$

The jacobian matrix for this system is given by:
$$J = \begin{pmatrix}
\frac{\partial \dot{R}}{\partial R} & \frac{\partial \dot{R}}{\partial E} \\ 
\frac{\partial \dot{E}}{\partial R} & \frac{\partial \dot{E}}{\partial E} 
\end{pmatrix}
=
\begin{bmatrix}
\alpha - \beta E & -\beta R \\
\delta E & -\gamma + \delta R
\end{bmatrix}$$

The stability of the fixed points can be determined by evaluating the eigenvalues of the jacobian matrix at each fixed point.

The eigenvalues at the first fixed point $(0, 0)$ are:

$$det(J-\lambda I) = 0 \implies (\alpha - \lambda)(-\gamma - \lambda) = 0$$
$$\implies \lambda_1 = \alpha = 2 , \quad \lambda_2 = -\gamma = -1 $$

This indicates that the first fixed point is a saddle point, which is unstable.

The eigenvalues at the second fixed point $\left(\frac{\gamma}{\delta}, \frac{\alpha}{\beta}\right)$ are:

$$det(J-\lambda I) = 0 \implies (\alpha - \beta \frac{\alpha}{\beta} - \lambda)(-\gamma + \delta \frac{\gamma}{\delta} - \lambda) - (-\beta \frac{\gamma}{\delta})(\delta \frac{\alpha}{\beta}) = 0$$
$$\implies (-\lambda)(-\lambda) - (-\beta \frac{\gamma}{\delta})(\delta \frac{\alpha}{\beta}) = 0$$
$$\implies \lambda^2 + \alpha \gamma = 0$$
$$\implies \lambda = \pm i \sqrt{\alpha \gamma} = \pm i \sqrt{2}$$

Since the eigenvalues are purely imaginary, the second fixed point is a center, which is stable but not asymptotically stable. This means that the system will exhibit oscillatory behavior around this fixed point.
