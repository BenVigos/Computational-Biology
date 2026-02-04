from scipy.integrate import solve_ivp 
from functools import partial
import matplotlib.pyplot as plt

# Vmax = 5.5 # mmol L-1 min-1
# Km = 50 # mM
# Ki = 3.5 # mM

Mlactose = 342.3 # g mol-1
milk_lactose_concentration = 25 # g L-1
volume_milk = 1000 # l

"""
First we want to convert all units to base units so mol, M (molar), g, liters, and minutes.
Vmax: 5.5 * 10^-3 / 60 = 5.5e-3 # mol L-1 min-1
Km: 50 * 10^-3 = 5.0e-2 # M
Ki: 3.5 * 10^-3 = 3.5e-3 # M
"""

Vmax = 5.5e-3 # mol L-1 min-1
Km = 5.0e-2 # M
Ki = 3.5e-3 # M

"""
v = vmax*S / Km + S = [mol L-1 min-1][mol] / [mol L-1] + [mol L-1] = [mol^2 L-1 min-1]/[2 (mol L-1)]
v = mol/min

dS/dt = v = (vmax * S)/(Km + s)
"""

safe_limit_concentration = 1 # g L-1

def lactase_activity(t, S, Vmax, Km):
    return -(Vmax * S)/(Km + S)

def inhibited_lactase_activity(t, y, Vmax, Km, Ki):
    S,I = y[0], y[1]
    return -(Vmax * S)/(Km * (1 + S/Ki) + S)

lactose_in_milk_M =  milk_lactose_concentration / Mlactose
safe_limit_M = safe_limit_concentration / Mlactose

def lactose_safe(t, S, Vmax, Km):
    return S[0] - safe_limit_M

lactose_safe.terminal = True
lactose_safe.direction = -1

### Without inhibition
print("Without inhibition:")

print(lactose_in_milk_M, safe_limit_M)
sol = solve_ivp(lactase_activity, (0, 100), [lactose_in_milk_M], args=(Vmax, Km), events=lactose_safe)

print(sol.y)
print(sol.y_events)
print(sol)

print(f"It took {sol.t_events[0][0]} minutes to reach a concentration below 1 g/L")
plt.plot(sol.t, sol.y[0])
plt.ylabel("Lactose concentration [mol]")
plt.xlabel("Time [min]")
plt.show()


### With inhibition
print("With inhibition:")