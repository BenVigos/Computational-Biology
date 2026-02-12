# Notes

## Equations and Model
Covers: Derivation of $K_{m2}$ for Question 1, the use of the Lineweaver-Burke plot (from which you can get Km1 and Vmax) for different S2.

Stoichiometric matrix N, steady state equillibrium.

Discus solver used (RK54), set max_step to 0.1. 

Key assumptions:

- Constant S2 in excess, because of that we use Km1 and Vmax values from S2=10000.
- Stoichiometric matrix describes reactants behaviour under different reactions v.

## Results and Figures

Lineweaver-Burke plot, table with Vmax and Km1 we read from that.

Km2 calculated by equation.
Km2 = 0.1000 +- 1.122e-11

Linear decrease of concentration over time.
"It took 728.58 seconds (12.1 minutes) to reach a concentration below 1 g/L."

Discuss the solution space plot of v1 vs v6

## Discussion

Lineweaver-Burke plot: since all lines are parrelel this is a ping pong mechanism. Since they all have the same slope and R^2 of 1 we can say this confidently.

It took 660.29 seconds (11.0 minutes) to reach a concentration below 1 g/L.

Linear line, expected exponential decay. Could be that due to high saturation of S2 we have always v of vmax which means the decay is linear with a rate vmax. This is likely also influenced by the design choice to keep S2 constant. We had to manually set max_step to 0.1 to prevent solve_ivp from going to negative concentrations.

Discuss $v_{6,max}\geq D$ and $v_{6,max}<D$ and changes in $V_\text{max}$ and $K_m$ affect feasible region.