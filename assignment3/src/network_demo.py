import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Network import NeuralNetwork, Pos

# --- 1. Setup & Run Simulation ---
grid_size = 10
u_grid = np.zeros((grid_size, grid_size))
u_grid[3:7, 3:7] = 1.5 
nn = NeuralNetwork(u_grid)

# Add neurons and connect
positions = [Pos(np.random.randint(0, 10), np.random.randint(0, 10)) for _ in range(15)]
for p in positions: nn.add_node(p.x, p.y)
for p1 in positions:
    for p2 in positions:
        if p1 != p2 and np.random.rand() < 0.4:
            nn.connect_nodes(p1, p2, w=0.5)

# Run for 100 steps and record weight snapshots
steps = 100
weight_history = []

for t in range(steps):
    if t % 20 == 0:
        for i in range(3): nn.G.nodes[positions[i]]["S"][-1] = 1
    
    V = nn.update_voltage(theta_e=0.05)
    S = nn.update_firing(V, theta_a=0.7, hill_n=6)
    nn.update_plasticity(eta=0.4, theta=0.6) # Higher eta for visible animation
    
    # Snapshot weights for animation
    current_weights = {(u, v): d['W'] for u, v, d in nn.G.edges(data=True)}
    weight_history.append(current_weights)

# --- 2. Animation Setup ---
fig, ax = plt.subplots(figsize=(8, 8))

def update(frame):
    ax.clear()
    ax.set_title(f"Neural Spatio-Temporal Plasticity: Step {frame}")
    
    # Draw U-Grid (Background)
    ax.imshow(u_grid.T, origin='lower', cmap='YlGn', alpha=0.3, extent=(-0.5, 9.5, -0.5, 9.5))
    
    # Draw Edges (Weights)
    weights = weight_history[frame]
    max_w = max(weights.values()) if weights else 1
    for (u, v), w in weights.items():
        linewidth = (w / max_w) * 2 if w > 0 else 0.1
        color = 'blue' if w > 0.5 else 'gray'
        ax.annotate("", xy=(v.x, v.y), xytext=(u.x, u.y),
                    arrowprops=dict(arrowstyle="->", color=color, alpha=0.3, lw=linewidth))

    # Draw Nodes (Firing State)
    for i, pos in enumerate(positions):
        is_firing = nn.G.nodes[pos]["S"][frame]
        color = 'red' if is_firing else 'black'
        size = 100 if is_firing else 30
        ax.scatter(pos.x, pos.y, c=color, s=size, zorder=5)
        ax.text(pos.x+0.1, pos.y+0.1, f"n{i}", fontsize=8)

    ax.set_xlim(-0.5, 9.5)
    ax.set_ylim(-0.5, 9.5)
    ax.set_xticks(range(10))
    ax.set_yticks(range(10))
    ax.grid(True, alpha=0.2)

ani = FuncAnimation(fig, update, frames=steps, interval=100, repeat=False)
plt.show()