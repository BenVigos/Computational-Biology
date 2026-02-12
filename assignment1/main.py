import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

kinetics = pd.read_csv("Kinetics.csv")

print('Question 1: The "Living" Ice')
print('Analyze the data and determine which mechanism "Europase" follows, justify yourconclusions statistically.')

print("|$S_2$|$K_m/V_\\text{max}$|$1/V\\text{max}$|$-1/K_{m1}$|$R^2$|$V\\text{max}$|$K_{m1}$|")
print("|-----|-------------------|----------------|-----------|-----|--------------|-----|")
for value in kinetics["S2"].unique():
    data = kinetics.loc[kinetics["S2"]==value]

    x = 1/data["S1"]
    y = 1/data["Rate"]
    n = len(data)

    plt.scatter(x, y, label=value)

    x_mean = x.mean()
    y_mean = y.mean()
    Sxy = np.sum(x*y) - n * x_mean * y_mean
    Sxx = np.sum(x*x) - n * x_mean * x_mean
    a = Sxy/Sxx
    b = y_mean-a*x_mean

    y_pred = a * x + b

    error = y - y_pred
    se = np.sum(error**2)
    mse = se/n 
    rmse = np.sqrt(mse)
    SSt = np.sum((y - y_mean)**2)
    R2 = 1- (se/SSt)

    print(f"|{value}|{a:.4f}|{b:.4f}|{-b/a:.4f}|{R2:.4f}|{1/b:.4f}|{-1/(-b/a):.4f}|")

    plt.axhline(y=0, color='k', linewidth=.5, alpha=0.5)
    plt.axvline(x=0, color='k', linewidth=.5, alpha=0.5)

    x_plot = np.arange(-b/a, max(x))
    y_plot = a * x_plot + b
    plt.plot(x_plot, y_plot)

plt.ylabel(r"$1/v$ (mM/s)")
plt.xlabel(r"$1/S_1$ (mM)")
plt.legend(title=r"$S_2$")
plt.title("Lineweaver-Burk plot")
plt.savefig("pingpong.jpg")
# plt.show()
plt.clf()

print("\nFind Km2")
S2 = 10000.0
Vmax = 1.0000
Km1 = 0.1000

Km2 = []
data = kinetics.loc[kinetics["S2"]==S2]
for index, row in data.iterrows():
    result = Vmax * S2 / row["Rate"] - Km1 * S2 / row["S1"] - S2
    Km2.append(result)

print(f"Km2 = {np.mean(Km2):.4f} +- {np.std(Km2):.3e}")

print()

# exit()
print("To synthesize the bio-polymer, an industrial batch reactor will be loaded with the co-factor (S2) in excess and an initial concentration of 100 g/L of the mineral substrate (S1). Calculate how much time (in seconds or minutes) is required to deplete the substrate S1 to below 1 g/L. Assume the Molecular Weight of substrate S1 is 150 g/mol.")

Vmax, Km1 = 1.0000, 0.100 # use values for [S2] = 10000.0 since S2 is in excess
S2 = 10000.0
Km2 = 0.1000
S1_MW = 150 # g/mol
S1_initial = 100 # g/L
S1_safe = 1 # g/L

S1_initial_mM = S1_initial / S1_MW * 1000
S1_safe_mM = S1_safe / S1_MW * 1000

def safe_limit(t, S, *args):
    return S[0] - S1_safe_mM

safe_limit.terminal = False
safe_limit.direction = -1

def reaction_ODE(t, S, Vmax, Km):
    """
    dS/dt = v = -(Vmax * S)/(Km + S)
    """

    return -(Vmax * S * S2) / (Km1*S2 + Km2*S + S*S2)


t_span = (0, 800)

sol = solve_ivp(
    reaction_ODE,
    t_span,
    [S1_initial_mM],
    args=(Vmax, Km1),
    events=[safe_limit],
    max_step=0.1,
)


print(f"It took {sol.t_events[0][0]:.2f} seconds ({sol.t_events[0][0]/60:.1f} minutes) to reach a concentration below {S1_safe} g/L.")
plt.title(r"Concentration of $S_1$ over time")
plt.plot(sol.t, sol.y[0])
plt.ylabel(r"$S_1$ (mM)")
plt.xlabel("Time (s)")
plt.axhline(y=S1_safe_mM, linestyle="dotted", alpha=0.5, color="k")
plt.plot(sol.t_events[0], sol.y_events[0][0], "x", color="k")

plt.savefig("Concentration.jpg")
# plt.show()
plt.clf()

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(sol.t, sol.y[0], label="Main Data")
ax.axhline(y=S1_safe_mM, linestyle="dotted", alpha=0.5, color="k")

if sol.t_events and len(sol.t_events[0]) > 0:
    ax.plot(sol.t_events[0], sol.y_events[0][0], "x", color="k", markersize=10)

ax.set_title(r"Concentration of $S_1$ over time")
ax.set_ylabel(r"$S_1$ (mM)")
ax.set_xlabel("Time (s)")

axins = ax.inset_axes([0.5, 0.5, 0.4, 0.4]) 

axins.plot(sol.t, sol.y[0])
axins.axhline(y=S1_safe_mM, linestyle="dotted", alpha=0.5, color="k")

x1, x2 = 665, 670
axins.set_xlim(x1, x2)

y_data_slice = sol.y[0][(sol.t >= x1) & (sol.t <= x2)]
if len(y_data_slice) > 0:
    y_min, y_max = y_data_slice.min(), y_data_slice.max()
    y_margin = (y_max - y_min) * 0.1
    axins.set_ylim(y_min - y_margin, y_max + y_margin)

axins.tick_params(labelsize=8)
# mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

plt.savefig("Concentration_with_Zoom.jpg", dpi=300)
