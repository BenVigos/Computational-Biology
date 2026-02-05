import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

kinetics = pd.read_csv("Kinetics.csv")

print('Analyze the data and determine which mechanism "Europase" follows, justify yourconclusions statistically.')

print("|$S_2$|Slope ($K_m/V_\\text{max}$)|Y-Intercept ($1/V\\text{max}$)|X-Intercept ($-1/K_m$)|$R^2$|$V\\text{max}$|$K_m$|")
print("|-----|---------------------------|------------------------------|----------------------|-----|--------------|-----|")
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
plt.title("Two substrate kinetics - Type 2 - Ping-pong mechanism")
plt.savefig("pingpong.jpg")
# plt.show()

print()

print("To synthesize the bio-polymer, an industrial batch reactor will be loaded with the co-factor (S2) in excess and an initial concentration of 100 g/L of the mineral substrate (S1). Calculate how much time (in seconds or minutes) is required to deplete the substrate S1 to below 1 g/L. Assume the Molecular Weight of substrate S1 is 150 g/mol.")

Vmax, Km = 1.0000, 0.1000 # use values for [S2] = 10000.0 since S2 is in excess
