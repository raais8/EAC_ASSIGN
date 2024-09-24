from sympy import Symbol, lambdify
from sympy.abc import x
from SymSolver import SymSolver
import numpy as np
import matplotlib.pyplot as plt

cParams = {
    'F': Symbol('F', real = True),
    'g': Symbol('g', real = True),
    'L': Symbol('L', positive = True),
    'A': Symbol('A', positive = True),
    'E': Symbol('E', positive = True),
    'r': Symbol('g', positive = True),  
}

SolverProb1 = SymSolver(cParams)
solution = SolverProb1.numSolve().simplify().evalf()
solution = solution.rhs
print(solution)

num_solution = lambdify(x, solution, 'numpy')
x_vals = np.linspace(0, 1, 1000)
y_vals = num_solution(x_vals)
plt.plot(x_vals, y_vals)
plt.show()
