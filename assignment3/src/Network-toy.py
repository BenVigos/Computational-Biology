import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

# --- 1. Simulation Parameters ---
GRID_SIZE = 15          
NEURON_DENSITY = 0.6    # Restored: 60% chance of a neuron existing
CONN_RADIUS = 3.5       # Restored: Max reach of connections
CONN_PROB = 0.3         # Restored: Chance to connect within radius
STEPS = 100             

THETA_E = 0.5           
HILL_N = 4.0            
ETA = 0.4               
THETA_U = 0.5           
W_MAX = 4.0             

# --- 2. Initialize the Sparse Grid Lattice ---
G = nx.DiGraph()

# Step A: Sparsely populate the grid
for x in range(GRID_SIZE):
    for y in range(GRID_SIZE):
        if np.random.rand() < NEURON_DENSITY:
            G.add_node((x, y), S=0, U=0.0, V=0.0)

# Step B: Create distance-based, sparse connections
nodes = list(G.nodes())
for i in range(len(nodes)):
    for j in range(len(nodes)):
        if i != j:
            n1, n2 = nodes[i], nodes[j]
            dist = np.sqrt((n1[0] - n2[0])**2 + (n1[1] - n2[1])**2)
            
            if dist <= CONN_RADIUS and np.random.rand() < CONN_PROB:
                G.add_edge(n1, n2, W=np.random.uniform(0.1, 0.2))

# Step C: Smart Input Node Selection
center_x, center_y = GRID_SIZE / 2, GRID_SIZE / 2
# Sort all nodes by their distance to the center
sorted_nodes = sorted(G.nodes(), key=lambda n: (n[0]-center_x)**2 + (n[1]-center_y)**2)

INPUT_NODE = None
# Find the closest node to the center that ACTUALLY has outgoing edges
for node in sorted_nodes:
    if G.out_degree(node) > 0:
        INPUT_NODE = node
        break

# Fallback just in case the random generation was incredibly sparse
if INPUT_NODE is None:
    INPUT_NODE = sorted_nodes[0]

pos = {node: node for node in G.nodes()} 

# --- 3. Diffusion Proxy ---
def simulate_fake_diffusion(graph):
    U_next = {}
    for node in graph.nodes():
        u_current = graph.nodes[node]['U']
        production = 2.0 if graph.nodes[node]['S'] == 1 else 0.0
        
        # Bleeds into connected neighbors
        neighbors = list(graph.successors(node)) + list(graph.predecessors(node))
        neighbor_u = sum(graph.nodes[nb]['U'] for nb in neighbors) / len(neighbors) if neighbors else 0.0
            
        u_new = (u_current * 0.4) + (neighbor_u * 0.4) + production
        u_new = (u_new + np.random.uniform(-0.05, 0.1)) * 0.85 
        
        U_next[node] = max(0.0, min(u_new, 3.0))
    return U_next

# --- 4. Animation Setup ---
plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(10, 10))

def update(frame):
    ax.clear() 
    
    S_new, W_new = {}, {}
    
    # Step A: Voltage and Firing Probability 
    for node in G.nodes():
        v_xy = sum([G.edges[pred, node]['W'] * G.nodes[pred]['S'] for pred in G.predecessors(node)])
        G.nodes[node]['V'] = v_xy
        
        if v_xy > 0:
            p_fire = (v_xy**HILL_N) / (THETA_E**HILL_N + v_xy**HILL_N)
        else:
            p_fire = 0.0
            
        if node == INPUT_NODE:
            p_fire = 1.0 
            
        S_new[node] = 1 if np.random.rand() < p_fire else 0
        
    nx.set_node_attributes(G, S_new, 'S') 
            
    # Step B: U Molecule Spread
    U_new = simulate_fake_diffusion(G)
    nx.set_node_attributes(G, U_new, 'U')

    # Step C: Spatio-Temporal Plasticity Update 
    for i, j in G.edges():
        w_t = G.edges[i, j]['W']
        s_i = G.nodes[i]['S']  
        u_j = G.nodes[j]['U']  
        
        w_next = w_t + ETA * s_i * (u_j - THETA_U)
        W_new[(i, j)] = max(0.0, min(w_next, W_MAX))

    nx.set_edge_attributes(G, W_new, 'W')

    # --- Draw the Network ---
    u_sizes = [G.nodes[n]['U'] * 400 + 10 for n in G.nodes()]
    s_colors = ['#ff0055' if G.nodes[n]['S'] == 1 else '#222222' for n in G.nodes()]
    
    edges = G.edges()
    weights = [G.edges[e]['W'] for e in edges]
    edge_widths = [w * 2.0 + 0.5 for w in weights] 

    sys.stdout.write(f"\rFrame {frame+1}/{STEPS} | Max Weight: {max(weights):.2f}")
    sys.stdout.flush()

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=s_colors, node_size=u_sizes, alpha=0.9, edgecolors='white')
    nx.draw_networkx_edges(
        G, pos, ax=ax, edgelist=edges, 
        width=edge_widths, 
        edge_color=weights, 
        edge_cmap=plt.cm.cool, 
        edge_vmin=0.0, edge_vmax=W_MAX, 
        alpha=0.6, arrows=False 
    )

    ax.set_title(f"Sparse Neural Network - Step {frame + 1}/{STEPS}", color='white', fontsize=14)
    ax.plot(INPUT_NODE[0], INPUT_NODE[1], 'wo', markersize=8) 
    ax.axis('equal')
    ax.axis('off')
    
ani = FuncAnimation(fig, update, frames=STEPS, interval=250, repeat=False)

plt.show(block=False)
plt.pause((STEPS * 250) / 1000.0 + 1.0)
plt.close()

print("\nSimulation complete.")