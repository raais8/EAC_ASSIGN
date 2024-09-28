from sympy import Symbol, lambdify
from sympy.abc import x
from SymSolver import ExactSolver
import numpy as np
import matplotlib.pyplot as plt

cParams = {
    'F': Symbol('F', real = True),
    'g': Symbol('g', real = True),
    'L': Symbol('L', positive = True),
    'A': Symbol('A', positive = True),
    'E': Symbol('E', positive = True),
    'r': Symbol('r', positive = True),  
}

SolverProb1 = ExactSolver(cParams)
solution = SolverProb1.numSolve().simplify().evalf()
print(solution)

num_solution = lambdify(x, solution, 'numpy')
x_vals = np.linspace(0, 1, 1000)
y_vals = num_solution(x_vals)
plt.plot(x_vals, y_vals,linewidth = 1.5)
plt.title(r'$\mathbf{Axial\, displacements}$', loc = 'center', fontsize=14)
plt.xlabel(r'$\mathbf{Position\,(m)}$', loc = 'center', fontsize = 12)
plt.ylabel(r'$\mathbf{Displacement\,(m)}$', loc = 'center', fontsize = 12)
plt.grid(True, linewidth = 1.5)
plt.tick_params(axis = 'both', which = 'major', labelsize = 11)
plt.gca().spines['top'].set_linewidth(1.5)
plt.gca().spines['right'].set_linewidth(1.5)  
plt.gca().spines['left'].set_linewidth(1.5) 
plt.gca().spines['bottom'].set_linewidth(1.5)
plt.xlim(0,1)
plt.ylim(-0.015,0.025)

plt.show()






















