import math
from sympy import Symbol, pi, Function, diff, dsolve, integrate
from sympy.abc import x
import numpy as np

class ExactSolver:

    def __init__(self, cParams):
        self.L, self.A, self.E, self.F = cParams['L'], cParams['A'], cParams['E'], cParams['F']
        self.g, self.r = cParams['g'], cParams['r']
        self.N = cParams['N']
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
        return(numSolution.rhs)
    
    def galerkinMethod(self):
        N_L = self.N.subs(x, self.L)
        B = diff(self.N, x)
        K1 = integrate(B.T * B, (x, 0, self.L))
        K2 = integrate((self.rho * self.N.T * self.N), (x, 0, self.L))
        K = K1 - K2
        F1 = integrate(self.N.T * self.s * self.r, (x, 0, self.L))
        F = pi**2 * self.g / self.L * N_L.T - F1
        r = [0]
        l = np.arange(1, self.N.shape[1]).tolist()
        KK = K.extract(l, l)
        FF = (F.extract(l,[0]) - K.extract(l, r) * self.g)
        dl = KK.LUsolve(FF)
        u = self.g + self.N.extract(l,[0]) * dl
        return(u)