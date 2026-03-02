class Cell:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.is_wall = False
        self.is_start = False
        self.is_goal = False
        self.g_cost = float('inf')  # Cost from start to this cell
        self.h_cost = float('inf')  # Heuristic cost to goal
        self.f_cost = float('inf')  # Total cost (g_cost + h_cost)
        self.parent = None  # Parent cell for path reconstruction

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        return hash((self.x, self.y))

    def __repr__(self):
        return f"Cell({self.x}, {self.y})"
