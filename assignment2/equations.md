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
du_i &= (c_i - \beta_i u_i)dt + \sigma_{1i}dW_1 \\
ds_i &= (b_i u_i - \gamma_i s_i)dt +\sigma_{2i}dW_2 \\
dp_i &= (k_{Pi} s_i - \delta_{Pi} p_i)dt
\end{align}$$


$$\begin{align}
du_i &= (c_i - \beta_i u_i)dt + \sigma_{1i}dW_1 \\
ds_i &= (b_i u_i - \gamma_i s_i)dt +\sigma_{2i}dW_2 \\
dp_i &= (k_{Pi} s_i - \delta_{Pi} p_i)dt
\end{align}$$


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

The jacobian matrix for this system is given by:
$$J = \begin{bmatrix}
\alpha - \beta E & -\beta R \\
\delta E & -\gamma + \delta R
\end{bmatrix}$$

The stability of the fixed points can be determined by evaluating the eigenvalues of the jacobian matrix at each fixed point.

The eigenvalues at the first fixed point $(0, 0)$ are:
$$\lambda_1 = \alpha = 2$$
$$\lambda_2 = -\gamma = -1$$