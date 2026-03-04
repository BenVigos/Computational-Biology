import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Simulation Parameters ---
GRID_SIZE = 15          # Larger grid to accommodate sparsity
NEURON_DENSITY = 0.4    # 40% chance a grid location has a neuron
CONN_RADIUS = 3.5       # Maximum distance an axon/dendrite can reach
CONN_PROB = 0.25        # 25% chance of forming a synapse if within radius
STEPS = 50

# Firing & Plasticity Parameters (kept from previous)
THETA_E = 1.5         
HILL_N = 4.0          
ETA = 0.05            
THETA_U = 0.5         
D_DIFF = 0.1          
DECAY_K = 0.05        
PROD_ALPHA = 0.8      
W_MAX = 3.0           

# --- 2. Initialize the Sparse Network ---
G = nx.DiGraph() # Directed graph for one-way synaptic connections

# Step A: Sparsely populate the grid with nodes
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
            # Calculate Euclidean distance between the two neurons
            dist = np.sqrt((n1[0] - n2[0])**2 + (n1[1] - n2[1])**2)
            
            # Connect if within radius AND passes the probability check
            if dist <= CONN_RADIUS and np.random.rand() < CONN_PROB:
                G.add_edge(n1, n2, W=np.random.uniform(0.1, 0.5))

# Step C: Assign the Input Node (find the neuron closest to the center)
center_x, center_y = GRID_SIZE / 2, GRID_SIZE / 2
distances_to_center = {n: (n[0]-center_x)**2 + (n[1]-center_y)**2 for n in G.nodes()}
INPUT_NODE = min(distances_to_center, key=distances_to_center.get)

# --- 3. Simulation Loop ---
for t in range(STEPS):
    S_new = {}
    U_new = {}
    W_new = {}
    
    # Calculate Voltage and Firing Probability
    for node in G.nodes():
        v_xy = sum([G.edges[pred, node]['W'] * G.nodes[pred]['S'] for pred in G.predecessors(node)])
        G.nodes[node]['V'] = v_xy
        
        p_fire = (v_xy**HILL_N) / (THETA_E**HILL_N + v_xy**HILL_N) if v_xy > 0 else 0.0
            
        if node == INPUT_NODE:
            p_fire = 1.0 # Force input
            
        S_new[node] = 1 if np.random.rand() < p_fire else 0
        
    # Calculate Diffusion (U Concentration)
    for node in G.nodes():
        u_t = G.nodes[node]['U']
        
        # Undirected neighbors for physical diffusion of the molecule
        neighbors = set(list(G.successors(node)) + list(G.predecessors(node)))
        laplacian = sum([G.nodes[nb]['U'] - u_t for nb in neighbors])
        
        p_xt = PROD_ALPHA * G.nodes[node]['S'] 
        U_new[node] = max(0.0, u_t + (D_DIFF * laplacian) - (DECAY_K * u_t) + p_xt)

    # Update Synaptic Weights
    for i, j in G.edges():
        w_t = G.edges[i, j]['W']
        s_i = G.nodes[i]['S']  
        u_j = G.nodes[j]['U']  
        
        w_next = w_t + ETA * s_i * (u_j - THETA_U)
        W_new[(i, j)] = max(0.0, min(w_next, W_MAX))

    # Apply updates
    nx.set_node_attributes(G, S_new, 'S')
    nx.set_node_attributes(G, U_new, 'U')
    nx.set_edge_attributes(G, W_new, 'W')

# --- 4. Visualization ---
plt.figure(figsize=(10, 10))
pos = {node: node for node in G.nodes()} # Use exact (x,y) grid coordinates

u_sizes = [G.nodes[n]['U'] * 150 + 20 for n in G.nodes()]
s_colors = ['#ff4d4d' if G.nodes[n]['S'] == 1 else '#a8d5e2' for n in G.nodes()]
edge_weights = [G.edges[e]['W'] * 1.5 for e in G.edges()]

# Draw the graph
nx.draw_networkx_nodes(G, pos, node_color=s_colors, node_size=u_sizes, alpha=0.9, edgecolors='black')
nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color='gray', alpha=0.4, arrows=True, arrowsize=10)

plt.title("Sparse Spatiotemporal Neural Network\nRed = Firing, Size = Plasticity ($U$), Line Thickness = Weight ($W$)")
plt.plot(INPUT_NODE[0], INPUT_NODE[1], 'go', markersize=12, label='Input Node (Center)')
plt.legend()
plt.axis('equal')
plt.axis('off')
plt.tight_layout()
plt.show()