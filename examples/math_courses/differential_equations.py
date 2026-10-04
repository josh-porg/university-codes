"""Differential equations homework 9 (``first semester/diff eq/diff_eq_homework_9``):
y'' + y' - 6y = 12 e^(-2t) and its homogeneous equation, solved with sympy."""

import sympy as sp

t = sp.Symbol("t")
y = sp.Function("y")
print(sp.dsolve(y(t).diff(t, 2) + y(t).diff(t) - 6 * y(t) - 12 * sp.exp(-2 * t)))
print(sp.dsolve(y(t).diff(t, 2) + y(t).diff(t) - 6 * y(t)))
