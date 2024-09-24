import math
from sympy import Symbol, pi, Function

class SymSolver:

    def __init__(self, r):
        self.x = Symbol('x')
        self.L = Symbol('L', positive=True)
        self.F = Symbol('F')
        self.E = Symbol('E', positive=True)
        self.A = Symbol('A', positive=True)
        self.g = Symbol('g')
        self.rho = pi / self.L**2
        self.s = self.g * pi**4 / self.L**2
        u = Function('u', real=True)
        self.q = self.E * (self.rho * u )



solver = SymSolver(2)