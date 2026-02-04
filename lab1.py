from scipy.integrate import solve_ivp 
from functools import partial
import matplotlib.pyplot as plt
import numpy as np

Vmax = 5.5 # mmol L-1 min-1
Km = 50 # mM
Ki = 3.5 # mM

Mlactose = 342.3 # g mol-1
milk_lactose_concentration = 25 # g L-1
volume_milk = 1000 # l

safe_limit_concentration = 1 # g L-1

def lactase_activity(t, S, Vmax, Km):
    """
    dS/dt = v = -(Vmax * S)/(Km + S)
    """
    return -(Vmax * S)/(Km + S)

lactose_in_milk_mM =  milk_lactose_concentration / Mlactose * 1000
safe_limit_mM = safe_limit_concentration / Mlactose * 1000

def lactose_safe(t, S, *args):
    return S[0] - safe_limit_mM

lactose_safe.terminal = True
lactose_safe.direction = -1

t_span = (0, 100)
t_eval = np.linspace(t_span[0], t_span[1], 1000)
### Without inhibition
print("Without inhibition:")

sol = solve_ivp(lactase_activity, t_span, [lactose_in_milk_mM], args=(Vmax, Km), events=lactose_safe, t_eval=t_eval)

print(f"It took {sol.t_events[0][0]} minutes to reach a concentration below 1 g/L")
plt.plot(sol.t, sol.y[0])
plt.ylabel("Lactose concentration [mol]")
plt.xlabel("Time [min]")
plt.axhline(y=safe_limit_mM, linestyle="dotted", alpha=0.5, color="k")
plt.show()


def lactase_activity_inhibited(t, y, Vmax, Km, Ki, S_initial):
    S = y[0]

    I = S_initial - S

    return -(Vmax*S)/(Km*(1+I/Ki)+S)



### With inhibition
print("With inhibition:")

t_span = (0, 500)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

sol = solve_ivp(lactase_activity_inhibited, t_span, [lactose_in_milk_mM], args=(Vmax, Km, Ki, lactose_in_milk_mM), events=lactose_safe, t_eval=t_eval)

print(f"It took {sol.t_events[0][0]} minutes to reach a concentration below 1 g/L")
plt.plot(sol.t, sol.y[0])
plt.ylabel("Lactose concentration [mol]")
plt.xlabel("Time [min]")
plt.axhline(y=safe_limit_mM, linestyle="dotted", alpha=0.5, color="k")
plt.show()
