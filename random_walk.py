import numpy as np

x_0 = (0, 0)
t = 7
N = 100

x_final = []

for _ in range(N):
    position = x_0
    for step in range(t):
        directions = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        direction = directions[np.random.choice(len(directions))]
        position = tuple(p + d for p, d in zip(position, direction))
    x_final.append(position)

x_final_arr = np.array(x_final)
x_0_arr = np.array(x_0)


squared_displacements = np.sum((x_final_arr - x_0_arr)**2, axis=1)
MSD = np.mean(squared_displacements)

""" 
MSD(t) = 2dDt
D = MSD(t)/2dt = MSD/(2*2*t)
"""

D = MSD/(2*len(x_0*t))

print(f"MSD: {MSD}; D: {D}")