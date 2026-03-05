# Import necessary FiPy modules
from fipy import CellVariable, Grid2D, TransientTerm, DiffusionTerm, Viewer, ImplicitSourceTerm

# --- 1. Simulation Parameters (The "Physics" and "Geometry") ---
nx = 50          # Number of grid points in x-direction
ny = 50          # Number of grid points in y-direction
dx = 1.          # Distance between grid points
dy = 1.
D = 1.0          # Nutrient diffusion coefficient
k = 0.5          # Nutrient consumption rate by the cell

# --- 2. The Domain (Our "Petri Dish") ---
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

# --- 3. The Variables (What we are modeling) ---
# Nutrient concentration variable
nutrient = CellVariable(name="Nutrient Concentration", mesh=mesh, value=0.0)

# A variable to mark where our cells are. 1.0 for a cell, 0.0 for empty space.
cell_locations = CellVariable(name="Cell Locations", mesh=mesh, value=0.0)

# --- 4. Boundary and Initial Conditions (The "Experiment Setup") ---

# A constant source of nutrient from the left wall (Dirichlet boundary condition)
# Concentration is fixed at 1.0 on the faces of cells at x=0.
nutrient.constrain(1.0, where=mesh.facesLeft)

# Place one cell in the middle of the dish
# Note: FiPy meshes are 1D arrays. To find the index from (x, y) coordinates:
# index = y * nx + x
center_index = int((ny / 2) * nx + (nx / 2))
cell_locations[center_index] = 1.0

# The other walls are impermeable (Neumann boundary condition)
# The flux across these faces is zero. This is the default, but we can be explicit.
# nutrient.faceGrad.constrain([0.], where=mesh.facesRight)
# nutrient.faceGrad.constrain([0.], where=mesh.facesTop)
# nutrient.faceGrad.constrain([0.], where=mesh.facesBottom)

# --- 5. The Governing Equation ---
# Consumption is proportional to the local nutrient concentration, there is also ImplicitSourceTerm to handle stiff equations
consumption_term = k * cell_locations * nutrient
eq = TransientTerm() == DiffusionTerm(coeff=D) - consumption_term

# --- 6. The Simulation Loop ---
# Setup the viewer to see the results
viewer = Viewer(vars=nutrient, datamin=0., datamax=1.)
viewer.plot()

timeStep = 0.2 # Ennsure CFL compliance
steps = 200

# Loop through time and solve the equation at each step
for step in range(steps):
    eq.solve(var=nutrient, dt=timeStep)
    viewer.plot()

print("Simulation finished. Press any key to close the viewer.")
input()