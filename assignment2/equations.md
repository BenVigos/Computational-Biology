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

## ODEs
$$\begin{align}
\frac{du_i}{dt} &= c_i - \beta_i u_i \\
\frac{ds_i}{dt} &= b_i u_i - \gamma_i s_i \\
\frac{dp_i}{dt} &= k_{Pi} s_i - \delta_{Pi} p_i
\end{align}$$

Applied to Mechanism 1, the ODEs are as follows:
$$\begin{align}
\frac{du_A}{dt} &= c_A - \beta_A u_A \\
\frac{ds_A}{dt} &= b_A u_A - \gamma_A s_A \\
\frac{du_B}{dt} &= c_B - \beta_B u_B \\
\frac{ds_B}{dt} &= b_B u_B - \gamma_B s_B \\
\frac{dp_A}{dt} &= k_{PA} s_A - \delta_{PA} p_A \\
\frac{dp_B}{dt} &= k_{PB} s_B - \delta_{PB} p_B
\end{align}$$

## SDEVelo

$$\begin{align}
du_i &= (c_i - \beta_i u_i)dt + \sigma_{1i}dW_1 \\
ds_i &= (b_i u_i - \gamma_i s_i)dt +\sigma_{2i}dW_2 \\
dp_i &= (k_{Pi} s_i - \delta_{Pi} p_i)dt
\end{align}$$

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
