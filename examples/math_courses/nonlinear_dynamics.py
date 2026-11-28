"""Math 650 nonlinear dynamics (``HW_1``, ``HW_2``, ``Logistic_Map_Unstable_Fixed_Points``,
``itterative_map_plotter``).
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

x, r, a, b, t = sp.symbols("x r a b t")

# HW 1: equilibria of x' = (1 - x) s x^a - x (1 - s)(1 - x)^a with s = 1/2, a = 5
f = (1 - x) * sp.Rational(1, 2) * x**5 - x * sp.Rational(1, 2) * (1 - x) ** 5
print("HW 1: f(0.1), f(0.9) =", float(f.subs(x, 0.1)), float(f.subs(x, 0.9)), " real equilibria:",
      sorted(float(s) for s in sp.solve(f, x) if s.is_real))

# HW 2
f2 = 3 * x**2 - 6 * x - 2
print("HW 2: roots of 3x^2 - 6x - 2:", sp.solve(f2, x))
g = x * (1 - x) * (2 - x)
print("      d/dx x(1-x)(2-x) =", sp.expand(sp.diff(g, x)), "; at x = 2:", sp.diff(g, x).subs(x, 2))
xd = -a * x * sp.log(b * x)
print("      x' = -a x ln(b x): equilibria", sp.solve(xd, x), "; slope at x = 1/b:", sp.simplify(sp.diff(xd, x).subs(x, 1 / b)))
k, m, g0 = sp.symbols("k m g", positive=True)
print("      0 = -2 k m x^2 + g ->", sp.solve(-2 * k * m * x**2 + g0, x))

# Logistic map: periodic points of f^n(x) = x for f = r x (1 - x)
fmap = lambda z: r * z * (1 - z)  # noqa: E731
it = x
for n in range(1, 4):
    it = fmap(it)
    print(f"Logistic map period-{n} points:", [sp.simplify(s) for s in sp.solve(sp.factor(it - x), x)][:6])

# Bifurcation diagram of x_{k+1} = x_k + beta - x_k^2
betas = np.arange(0.001, 2.0, 0.002)
pb, px = [], []
for beta in betas:
    z = 0.5
    with np.errstate(over="ignore", invalid="ignore"):
        for _ in range(2000):
            z = z + beta - z**2
        z_ss = z
        for _ in range(500):
            z = z + beta - z**2
            pb.append(beta)
            px.append(z)
            if abs(z - z_ss) < 1e-3:
                break
fig, ax = plt.subplots(facecolor="k")
ax.plot(pb, px, ".", ms=0.8, color="w")
ax.set(facecolor="k", xlabel="beta", ylabel="x", ylim=(-1, 3), title="x + beta - x^2")
plt.show()
