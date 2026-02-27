---
geometry:
- margin=1in
---
# Assignment 2: Model development and equations 

## Viterbi algorithm
The Viterbi algorithm has the following steps:
1. *Initialization:* Set the initial probabilities for each state at time t=0.
2. *Recursion:* For each time step t and each state, calculate the maximum probability of being in that state given the previous states and the observed data. This involves using the transition probabilities and the emission probabilities.
3. *Termination:* After processing all time steps, identify the state with the highest probability at the final time step.
4. *Backtracking:* Trace back through the states to determine the most likely sequence of states that led to the observed data.

## Question 1

### ODEs

$$
\begin{aligned}
\frac{dr_a}{dt} &= m_a \times h^+(P_b, \theta_b, n_b) + \gamma_a \times r_a \\
\frac{dr_b}{dt} &= m_b \times h^-(P_a, \theta_a, n_a) + \gamma_b \times r_b \\
\frac{dP_a}{dt} &= k_{Pa} \times r_a - \delta_{Pa} \times P_a \\  
\frac{dP_b}{dt} &= k_{Pb} \times r_b - \delta_{Pb} \times P_b
\end{aligned}
$$

When transcriptional hijack, $h^-(P_a, \theta_a, n_a) = 1$

| Parameter | Definition |
| :--- | :--- |
| $dr_i$ | Change in concentration of transcribed mRNA A or B. |
| $m_i$ | Maximum transcription rate coefficient for mRNA ($Ms^{-1}$). |
| $h^+(P, \theta, n)$ | Hill activation function. |
| $h^-(P, \theta, n)$ | Hill inhibition function. |
| $\theta_i$ | Expression threshold for protein binding ($M$). |
| $n_i$ | Hill coefficient representing regulatory nonlinearity. |
| $\gamma_i$ | Degradation rate of mRNA ($s^{-1}$). |
| $r_i$ | Concentration of transcribed mRNA A or B ($M$). |
| $dP_i$ | Change in concentration of Protein A or Protein B. |
| $k_{Pi}$ | Translation rate of Protein A or Protein B ($s^{-1}$). |
| $P_i$ | Concentration of Protein ($M$). |
| $\delta_{Pi}$ | Degradation rate of Protein A or Protein B ($s^{-1}$). |
| $dt$ | Time differential. |

### SDEVelo

$$
\begin{aligned}
\alpha_A(t) &= \frac{c_A}{1+\exp b_A(t-a_A)} \\
\beta_A^* &= \beta_A h^+(P_B, \theta_B, n_B) \\
dU_A &= (\alpha_A(t) - \beta_A^* U_A(t))dt + \sigma_{1A}dW_{1A} \\
dS_A &= (\beta_A^* U_A(t) - \gamma_A S_A(t))dt + \sigma_{2A}dW_{2A} \\
dP_A &= (k_{PA} S_A(t) - \delta_{PA} P_A(t))dt \\
\end{aligned}
$$

$$
\begin{aligned}
\alpha_B(t) &= \frac{c_B}{1+\exp b_B(t-a_B)} \\
\beta_B^* &= \beta_B h^-(P_A, \theta_A, n_A) \\
dU_B &= (\alpha_B(t) - \beta_B^* U_B(t))dt + \sigma_{1B}dW_{1B} \\
dS_B &= (\beta_B^* U_B(t) - \gamma_B S_B(t))dt + \sigma_{2B}dW_{2B} \\
dP_B &= (k_{PB} S_B(t) - \delta_{PB} P_B(t))dt \\
\end{aligned}
$$

Table: Definitions of the parameters used in the SDEVelo equations

| Parameter | Definition |
| :--- | :--- |
| $\alpha_i(t)$ | Time-dependent transcription rate of unspliced mRNA (pre-mRNA). |
| $c_i$ | Maximum transcription rate coefficient for pre-mRNA ($Ms^{-1}$). |
| $b_i$ | Steepness parameter of the transcription rate sigmoid function. |
| $t$ | Time variable. |
| $a_i$ | Time shift or activation delay for transcription ($s$). |
| $\beta^*_i$ | Regulated splicing rate influenced by protein interactions. |
| $\beta_i$ | Base splicing rate parameter ($s^{-1}$). |
| $h^+(P, \theta, n)$ | Hill activation function. |
| $h^-(P, \theta, n)$ | Hill inhibition function. |
| $P_i$ | Concentration of Protein ($M$). |
| $\theta_i$ | Expression threshold for protein binding ($M$). |
| $n_i$ | Hill coefficient representing regulatory nonlinearity. |
| $dU_i$ | Change in concentration of unspliced mRNA (pre-mRNA). |
| $U_i(t)$ | Concentration of unspliced mRNA (pre-mRNA) at time $t$ ($M$). |
| $dt$ | Time differential. |
| $\sigma_{1i}$ | Noise intensity parameter for pre-mRNA transcription ($Ms^{-1/2}$). |
| $dW_{1i}$ | Differential of the Wiener process for transcription noise. |
| $dS_i$ | Change in concentration of spliced mRNA. |
| $S_i(t)$ | Concentration of spliced mRNA at time $t$ ($M$). |
| $\gamma_i$ | Degradation rate of spliced mRNA ($s^{-1}$). |
| $\sigma_{2i}$ | Noise intensity parameter for the splicing process ($Ms^{-1/2}$). |
| $dW_{2i}$ | Differential of the Wiener process for splicing noise. |
| $dP_i$ | Change in concentration of Protein A or Protein B. |
| $k_{Pi}$ | Translation rate of Protein A or Protein B ($s^{-1}$). |
| $\delta_{Pi}$ | Degradation rate of Protein A or Protein B ($s^{-1}$). |


## Downstream metabolic effects
$$
\begin{aligned}
\frac{dR}{dt} &= \alpha R - \beta RE \\
\frac{dE}{dt} &= -\gamma E + \delta RE
\end{aligned}
$$

where:
* $\alpha R$: the growth of the resource R, which is proportional to its current amount. This term represents the natural growth or replenishment of the resource in the absence of any interactions with the enzyme E.
* $-\beta RE$: the consumption of the resource R by the enzyme E. This term represents the rate at which the enzyme E utilizes the resource R, and it is proportional to both the amount of resource R and the amount of enzyme E.
* $−\gamma E$: the natural decay or degradation of the enzyme E. This term represents the rate at which the enzyme E is lost or deactivated over time, independent of its interaction with the resource R.
* $\delta RE$: the production or activation of the enzyme E by the resource R. This term represents the rate at which the enzyme E is generated or activated in response to the presence of the resource R, and it is proportional to both the amount of resource R and the amount of enzyme E.

with $\alpha=2$, $\beta=1.1$, $\gamma=1$, $delta=0.9$, $R(0)=1$ and $E(0)=0.5$

The fixed points of this system can be found by setting the derivatives to zero:
$$
\begin{aligned}
0 &= \alpha R - \beta R E \\
0 &= -\gamma E + \delta E R
\end{aligned}
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
$$\implies \lambda_1 = \alpha = 2, \quad \lambda_2 = -\gamma = -1$$

This indicates that the first fixed point is a saddle point, which is unstable.

The eigenvalues at the second fixed point $\left(\frac{\gamma}{\delta}, \frac{\alpha}{\beta}\right)$ are:

$$det(J-\lambda I) = 0 \implies (\alpha - \beta \frac{\alpha}{\beta} - \lambda)(-\gamma + \delta \frac{\gamma}{\delta} - \lambda) - (-\beta \frac{\gamma}{\delta})(\delta \frac{\alpha}{\beta}) = 0$$
$$\implies (-\lambda)(-\lambda) - (-\beta \frac{\gamma}{\delta})(\delta \frac{\alpha}{\beta}) = 0$$
$$\implies \lambda^2 + \alpha \gamma = 0$$
$$\implies \lambda = \pm i \sqrt{\alpha \gamma} = \pm i \sqrt{2}$$

Since the eigenvalues are purely imaginary, the second fixed point is a center, which is stable but not asymptotically stable. This means that the system will exhibit oscillatory behavior around this fixed point.
