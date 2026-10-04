"""AE 573 thermodynamics HW 18: property changes of a solid with the equation of state
``p = (a + b T) v / (1 + (((a + b T) eps0 / theta0 - 2) / eps0) v)`` (``Thermo_HW_18``, ``GasLawHW18``).

``v`` is the strain, ``p`` the stress. The specific heats follow the
homework: ``c_v = T (dv/dT)_p^2 (dp/dv)_T / (rho (k - 1))`` with
``c_p = k c_v``; the changes of h, u and s between 330 K and 360 K are
integrated at the fixed strain ``v = 1.75e-6``. ``GasLaw_Test`` (ideal-gas
and Redlich-Kwong checks of ``unicodes.thermo.GasLaw``) is covered by the tests.
"""

import sympy as sp
from scipy.integrate import quad

# The implicit derivatives are formed directly: building a full GasLaw (with its symbolic
# integrals) for this equation of state takes minutes.
T, v, p, a, b, eps0, th0 = sp.symbols("T v p a b epsilon_0 theta_0")
k, rho = sp.symbols("k rho")
F = p - (a + b * T) * v / (1 + (((a + b * T) * eps0 / th0 - 2) / eps0) * v)


def d(x, y):
    return -sp.diff(F, y) / sp.diff(F, x)


P = {"VrTcP": d(v, T), "PrVcT": d(p, v), "TrPcV": d(T, p), "TrVcP": d(T, v), "VrPcT": d(v, p), "PrTcV": d(p, T)}
values = {sp.Symbol("epsilon_0"): 2500e-6, sp.Symbol("theta_0"): 6e7, sp.Symbol("a"): 31.9e9, sp.Symbol("b"): -1e7,
          k: 1.1, rho: 1800, v: 1.75e-6}
p_of_T = sp.solve(F, p)[0]
cv = (T * P["VrTcP"] ** 2 * P["PrVcT"] / (rho * (k - 1))).subs(p, p_of_T)
cv_f = sp.lambdify(T, cv.subs(values))
cp_f = sp.lambdify(T, (k * cv).subs(values))
print("reciprocity (dT/dv)(dv/dT) = 1:", sp.simplify(P["TrVcP"] * P["VrTcP"] - 1) == 0)
print(f"stress at 300 K: {float(p_of_T.subs(values).subs(T, 300)):.4g}")
print(f"(dT/dp)_v at 300 K: {float(P['TrPcV'].subs(values).subs(p, p_of_T.subs(values)).subs(T, 300)):.4g}")
print(f"c_v(300 K) = {cv_f(300):.4g}, c_p(300 K) = {cp_f(300):.4g}")
dh = quad(cp_f, 330, 360)[0]
du = quad(cv_f, 330, 360)[0]
ds = quad(lambda t: cp_f(t) / t, 330, 360)[0]
print(f"330 -> 360 K: delta h = {dh:.4g}, delta u = {du:.4g}, delta s = {ds:.4g}")
