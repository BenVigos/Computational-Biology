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

$$\begin{align}
\frac{dr_a}{dt} &= m_a \times h^+(P_b, \theta_b, n_b) + \gamma_a \times r_a \\
\frac{dr_b}{dt} &= m_b \times h^-(P_a, \theta_a, n_a) + \gamma_b \times r_b \\
\frac{dP_a}{dt} &= k_a \times r_a - \delta_a \times P_a \\  
\frac{dP_b}{dt} &= k_b \times r_b - \delta_b \times P_b
\end{align}$$

where $h^+(P_b, \theta_b, n_b)$ and $h^-(P_a, \theta_a, n_a)$ are hill activation and inhibition function respectively. When transcriptional hijack, $h^-(P_a, \theta_a, n_a) = 1$

### SDEVelo

$$\begin{align}
\alpha_A(t) &= \frac{c_A}{1+\exp b_A(t-a_A)} \\
\beta_A^* &= \beta_A h^+(P_B, \theta_B, n_B) \\
dU_A &= (\alpha_A(t) - \beta_A^* U_A(t))dt + \sigma_{1A}dW_{1A} \\
dS_A &= (\beta_A^* U_A(t) - \gamma_A S_A(t))dt + \sigma_{2A}dW_{2A} \\
dP_A &= (k_{PA} S_A(t) - \delta_{PA} P_A(t))dt \\
\end{align}$$

$$\begin{align}
\alpha_B(t) &= \frac{c_B}{1+\exp b_B(t-a_B)} \\
\beta_B^* &= \beta_B h^-(P_A, \theta_A, n_A) \\
dU_B &= (\alpha_B(t) - \beta_B^* U_B(t))dt + \sigma_{1B}dW_{1B} \\
dS_B &= (\beta_B^* U_B(t) - \gamma_B S_B(t))dt + \sigma_{2B}dW_{2B} \\
dP_B &= (k_{PB} S_B(t) - \delta_{PB} P_B(t))dt \\
\end{align}$$

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

## Question 2

$$\begin{align}
\frac{dR}{dt} &= \alpha R - \beta RE \\
\frac{dE}{dt} &= -\gamma E + \delta RE\\
\end{align}$$

where:
* $\alpha R$:
* $-\beta RE$:
* $−\gamma E$:
* $\delta RE$:

We use  $\alpha = 2, \beta = 1.1, \gamma = 1$ and $\delta = 0.9$.

For stability analysis, first we find stationary points as: 
$\frac{dR}{dt} = R(\alpha - \beta E) = 0$ and
$\frac{dE}{dt} = E(-\gamma + \delta R) = 0$
obtaining:
* $R_1^* = 0$, $E_1^* = 0$
* $R_2^* = \frac{\gamma}{\delta} = \frac{1}{0.9} \approx 1.11$,  $E_2^* = \frac{\alpha}{\beta} = \frac{2}{1.1} \approx 1.81$

Then, Jacobian is:
$$ J = \begin{pmatrix} \frac{\partial \dot{R}}{\partial R} & \frac{\partial \dot{R}}{\partial E} \\ \frac{\partial \dot{E}}{\partial R} & \frac{\partial \dot{E}}{\partial E} \end{pmatrix} = \begin{pmatrix} \alpha - \beta E & -\beta R \\ \delta E & -\gamma + \delta R \end{pmatrix}$$

At point $R=0, E=0$: $$J(0,0) = \begin{pmatrix} 2 & 0 \\ 0 & -1 \end{pmatrix}$$

So that $\lambda_1 = 2$ e $\lambda_2 = -1$, then $(0,0)$ is point of saddle (unstable).

At point, $(R = 1.11, E = 1.81)$, $$J(R^*,E^*) = \begin{pmatrix} 0 & -1.22 \\ 1.63 & 0 \end{pmatrix}$$.

Using $\lambda^2 - \text{trace}(J)\lambda + \text{det}(J) = 0$, with $\text{det}(J) \approx 2$, then  $\lambda^2 + 2 = 0$, from which: $\lambda_{1,2} = \pm i \sqrt{2}$. Since eigenvalues are purely imaginary, the equilibrium is a center: stable, not asymptotically stable.

## Statistical measures



## Downstream metabolic effects
$$\begin{align}
\frac{dR}{dt} = \alpha R - \beta R E\\
\frac{dE}{dt} = -\gamma E + \delta E R
\end{align}$$

with $\alpha=2$, $\beta=1.1$, $\gamma=1$, $delta=0.9$, $R(0)=1$ and $E(0)=0.5$

The fixed points of this system can be found by setting the derivatives to zero:
$$\begin{align}
0 &= \alpha R - \beta R E \\
0 &= -\gamma E + \delta E R
\end{align}$$

and solving for R and E. This yields the fixed points:
1. $(R^*, E^*) = (0, 0)$
2. $(R^*, E^*) = \left(\frac{\gamma}{\delta}, \frac{\alpha}{\beta}\right)$
3. $(R^*, E^*) = \left(\frac{\alpha}{\beta}, \frac{\gamma}{\delta}\right)$

The jacobian matrix for this system is given by:
$$J = \begin{bmatrix}
\alpha - \beta E & -\beta R \\
\delta E & -\gamma + \delta R
\end{bmatrix}$$
