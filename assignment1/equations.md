---
geometry:
- margin=1in
---

# Assignment 1: Equations
##  Question 1: The "Living" Ice
### Derivation of $K_{m2}$
$$\begin{align}
v &= \frac{V_{max}S_1S_2}{K_{m1}S_2+K_{m2}S_1+S_1S_2} \\
v(K_{m1}S_2+K_{m2}S_1+S_1S_2) &= V_{max}S_1S_2 \\
K_{m1}S_2+K_{m2}S_1+S_1S_2 &= \frac{V_{max}S_1S_2}{v} \\
K_{m2}S_1 &= \frac{V_{max}S_1S_2}{v} - K_{m1}S_2 - S_1S_2 \\
K_{m2} &= \frac{V_{max}S_2}{v} - \frac{K_{m1}S_2}{S_1} - S_2
\end{align}$$

## Question 2: The Case of the possible Biomass
$$N=\begin{bmatrix}
-1&&&&&&&1&-1 \\
1&-1&&&&&&& \\
&1&-1&&&&&& \\
&&1&-1&&&&& \\
&&&1&-1&&&& \\
&&&&1&-1&&& \\
&&&&&1&-1&& \\
&&&&&&1&-1& \\
&&&&&&&&1 \\
\end{bmatrix}$$