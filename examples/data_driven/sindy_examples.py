"""SINDy on the Lorenz and Van der Pol systems (``SINDy_Version_0`` ... ``SINDy_Version_2_6``,
``simLorenz``, ``simVanderpol``).

Version 0 regresses exact Lorenz derivatives; versions 1-2.6 identify the
Van der Pol oscillator from total-variation regularised derivatives
(``TVRegDiff``), then simulate the identified model from a new initial
condition::

    python sindy_examples.py --system lorenz
    python sindy_examples.py --system vanderpol --tv
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from unicodes.decomposition import SINDy, lorenz, van_der_pol
from unicodes.numerics import tv_derivative

p = argparse.ArgumentParser()
p.add_argument("--system", choices=["lorenz", "vanderpol"], default="vanderpol")
p.add_argument("--tv", action="store_true", help="estimate derivatives with TV regularisation")
p.add_argument("--threshold", type=float, default=0.5)
p.add_argument("--ridge", type=float, default=0.0)
p.add_argument("--t-end", type=float, default=20.0)
a = p.parse_args()

dt = 0.001
t = np.arange(dt, a.t_end, dt)
if a.system == "lorenz":
    f, x0, names, terms, x0_test = lorenz, [0, 1, 20], ["x", "y", "z"], ("const", "linear", "poly2"), [1, 1, 25]
else:
    f, x0, names, terms, x0_test = van_der_pol, [0, 1], ["x", "y"], ("const", "linear", "poly2", "poly3"), [0.5, 0.5]
X = solve_ivp(f, (t[0], t[-1]), x0, t_eval=t, rtol=1e-12, atol=1e-12).y.T
dX_exact = np.array([f(0, x) for x in X])
if a.tv:
    # TVRegDiff(x, 10, 2e-5, [], "small", 1e12, dt): returns n + 1 edge values, the last is dropped
    dX = np.column_stack([tv_derivative(X[:, k], 10, 2e-5, dt, ep=1e12)[:-1] for k in range(X.shape[1])])
else:
    dX = dX_exact
model = SINDy(terms, threshold=a.threshold, ridge=a.ridge).fit(X, dX)
print("\n".join(model.equations(names)))

X_true = solve_ivp(f, (t[0], t[-1]), x0_test, t_eval=t, rtol=1e-12, atol=1e-12).y.T
X_model = model.simulate(x0_test, t)
fig = plt.figure()
ax = fig.add_subplot(projection="3d" if len(names) == 3 else None)
ax.plot(*X_true.T, "k", lw=1, label="true")
ax.plot(*X_model.T, "r--", lw=1, label="SINDy model")
ax.legend()
fig, ax = plt.subplots()
ax.plot(dX_exact[:, 0], dX_exact[:, 1], "y", label="exact derivative")
ax.plot(dX[:, 0], dX[:, 1], "k", lw=0.5, label="derivative used")
ax.plot(*model.predict(X)[:, :2].T, "m--", label="model")
ax.legend()
plt.show()
