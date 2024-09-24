
import math
from sympy import Symbol
from sympy.abc import x
from SymSolver import SymSolver

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
print(solution)