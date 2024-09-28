import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# ================================
# 1. Analytical Solution using SymPy
# ================================

# Define the variable and the function for the analytical solution
x = sp.symbols('x')
C1, C2 = sp.symbols('C1 C2')

# General solution to the ODE y''(x) - 2 = 0
y_general = -x**2 + C1*x + C2
print(f"General solution: {y_general}")

# Apply boundary condition y(0) = -1 to find C2
C2_value = sp.solve(y_general.subs(x, 0) + 1, C2)[0]

# Substitute C2 into the general solution
y_general = y_general.subs(C2, C2_value)

# Apply boundary condition y'(1) = 0 to find C1
y_prime = sp.diff(y_general, x)
C1_value = sp.solve(y_prime.subs(x, 1), C1)[0]

# Substitute C1 into the general solution
exact_solution = y_general.subs(C1, C1_value)
print(f"Exact solution: {exact_solution}")

# Convert the exact solution into a lambda function for numerical use
exact_sol_lambdified = sp.lambdify(x, exact_solution, 'numpy')

# ================================
# 2. FEM Solution with 4 nodes
# ================================

# Set the number of nodes and elements
N = 4
L = 1  # Domain length
nodes = np.linspace(0, L, N)
h = nodes[1] - nodes[0]  # Step size

# Global stiffness matrix K and load vector f
K = np.zeros((N, N))
f = np.zeros(N)

# Define the source term g(x) = -2
g = -2

# Assemble the stiffness matrix K and load vector f
for i in range(1, N - 1):
    K[i, i] += 2 / h
    K[i, i - 1] -= 1 / h
    K[i, i + 1] -= 1 / h
    f[i] += g * h  # Correct the load vector for the internal nodes

# Neumann boundary at x = 1: y'(1) = 0
K[N - 1, N - 1] += 1 / h
K[N - 1, N - 2] -= 1 / h
f[N - 1] += g * h / 2  # Apply half the load for Neumann node at x = 1

# Dirichlet boundary at x = 0: y(0) = -1
K[0, 0] = 1
K[0, 1] = 0  # Clear out any off-diagonal entries
f[0] = -1  # Apply the boundary condition directly in the load vector

# Solve the linear system
Y = np.linalg.solve(K, f)

# ================================
# 3. Plot the Exact and FEM Solutions
# ================================

# Plotting the exact solution
x_values = np.linspace(0, L, 100)
y_exact_values = exact_sol_lambdified(x_values)

# Plot the FEM solution using piecewise linear interpolation
y_fem_values = np.zeros_like(x_values)
for i in range(N - 1):
    # Interpolate between each pair of FEM nodes
    mask = (x_values >= nodes[i]) & (x_values <= nodes[i + 1])
    y_fem_values[mask] = np.interp(x_values[mask], nodes[i:i + 2], Y[i:i + 2])

# Create the plot
plt.plot(x_values, y_exact_values, label='Exact Solution', color='blue')
plt.plot(x_values, y_fem_values, label='FEM Solution', color='red', linestyle='--')

# Add labels, legend, and title
plt.xlabel('x')
plt.ylabel('y(x)')
plt.legend()
plt.title('Exact Solution vs FEM Solution')
plt.grid(True)

# Show the plot
plt.show()
