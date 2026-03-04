import networkx as nx
import numpy as np
from typing import NamedTuple


class Pos(NamedTuple):
    x: int
    y: int


class NeuralNetwork:
    def __init__(self, grid: np.ndarray):
        self.grid = grid  # 2D array for U_{x,y}
        self.G = nx.DiGraph()

    def add_node(self, x: int, y: int):
        pos = Pos(x, y)
        # Initialize history with starting values
        self.G.add_node(pos, S=[0], v=[0.0])

    def connect_nodes(self, n1_pos: Pos, n2_pos: Pos, w: float):
        """Adds a weighted edge between two existing nodes."""
        self.G.add_edge(n1_pos, n2_pos, W=w)

    def update_voltage(self, theta_e: float):
        """
        $$v_{x,y}(t)=\sum_{i,j}W_{(i,j)\rightarrow(x,y)}(t)\cdot(S_{i,j}(t)-\theta_e)$$
        """
        nodes = list(self.G.nodes())

        # Get firing states and subtract threshold: (S_j - theta_e)
        S_t = np.array([self.G.nodes[n]["S"][-1] for n in nodes], dtype=float)

        # Get weight matrix (W_ij)
        W = nx.to_numpy_array(self.G, nodelist=nodes, weight="W")

        # Matrix multiply: Rows are sources, columns are targets
        # We need S_shifted @ W to sum over sources for each target
        V_t = (S_t - theta_e) @ W

        for i, n in enumerate(nodes):
            self.G.nodes[n]["v"].append(V_t[i])

        return V_t

    def update_firing(self, V_t: np.ndarray, theta_a: float, hill_n: float):
        """
        $$P(S_{x,y}(t+1)=1)=\frac{v_{x,y}(t)^n}{\theta_a^n+v_{x,y}(t)^n}$$
        """
        # Rectify V to handle negative results from (S - theta_e)
        v_pos = np.maximum(V_t, 0)
        p_fire = (v_pos**hill_n) / (theta_a**hill_n + v_pos**hill_n)
        S_next = (np.random.rand(len(V_t)) < p_fire).astype(int)

        nodes = list(self.G.nodes())
        for i, n in enumerate(nodes):
            self.G.nodes[n]["S"].append(S_next[i])

        return S_next

    def update_plasticity(self, eta: float, theta: float):
        """
        $$W_{(i,j)\rightarrow(x,y)}(t+1)=W_{(i,j)\rightarrow(x,y)}(t)+\eta\cdot S_{i,j}(t)\cdot(U_{x,y}(t)-\theta)$$
        """
        for u, v, data in self.G.edges(data=True):
            # i,j = source (u) | x,y = target (v)
            S_ij = self.G.nodes[u]["S"][-1]
            U_xy = self.grid[v.x, v.y]

            data["W"] += eta * S_ij * (U_xy - theta)
