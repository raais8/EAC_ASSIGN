import math
from sympy import Symbol, pi, Function, diff, dsolve
from sympy.abc import x

class SymSolver:

    def __init__(self, cParams):
        self.L, self.A, self.E, self.F = cParams['L'], cParams['A'], cParams['E'], cParams['F']
        self.g, self.r = cParams['g'], cParams['r']
        self.u = Function('u', real=True)
        self.rho = pi**3 / self.L**2
        self.s = self.g * pi**4 / self.L**2
        self.q = self.E * (self.rho * self.u(x) - self.s * self.r)
    
    def buildAndSolveODE(self):
        up = diff(self.u(x),x)
        upp = diff(self.u(x), x, x)
        f = self.q / self.E
        EDO = upp + f
        EDO = EDO.subs(self.r, (x / self.L)**2)
        solution = dsolve(EDO, self.u(x), ics={self.u(0): -self.g, self.u(x).diff(x).subs(x,self.L): pi**2 * self.g / self.L})
        return(solution)
    
    def numSolve(self):
        solution = self.buildAndSolveODE()
        numSolution = solution.subs({self.L: 1, self.g: 0.01})
        return(numSolution)
