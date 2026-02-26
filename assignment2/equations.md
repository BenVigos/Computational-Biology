---
geometry:
- margin=1in
---
# Assignment 2: Model development and equations 

## ODEs

$$\begin{align}
\frac{dr_a}{dt} &= m_a * h^+(P_b, \theta_b, n_b) + \gamma_a * r_a \\
\frac{dr_b}{dt} &= m_b * h^-(P_a, \theta_a, n_a) + \gamma_b * r_b \\
\frac{dP_a}{dt} &= k_a * r_a - \delta_a * P_a \\  
\frac{dP_b}{dt} &= k_b * r_b - \delta_b * P_b
\end{align}$$

where h^+(P_b, \theta_b, n_b) and h^-(P_a, \theta_a, n_a) are hill activation and inhibition function respectively. When transcriptional hijack, $h^-(P_a, \theta_a, n_a) = 1$

## SDEVelo

$$\begin{align}
du_i &= (c_i - \beta_i u_i)dt + \sigma_{1i}dW_1 \\
ds_i &= (b_i u_i - \gamma_i s_i)dt +\sigma_{2i}dW_2 \\
dp_i &= (k_{Pi} s_i - \delta_{Pi} p_i)dt
\end{align}$$
