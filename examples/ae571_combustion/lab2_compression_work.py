"""AE 571 lab 2: compression and expansion work from a simulated adiabatic piston cycle
(``AE_571_lab_2_Poznanski``, ``AE_571_lab_2_Extra_credit_Poznanski``).

Input is the pair of Fluent monitor files
(``volume-averaged-density-rfile_6_1.out``, ``volume-averaged-pressure-rfile_6_1.out``)
or, for the extra credit, ``AE571_Lab02_Simulation.mat`` (``va_density``,
``va_pressure``, ``time``)::

    python lab2_compression_work.py --out-files density.out pressure.out --cp 662.66
    python lab2_compression_work.py --mat AE571_Lab02_Simulation.mat --cp 1200 --split 1000
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, trapezoid

from unicodes.io import load_mat

p = argparse.ArgumentParser()
g = p.add_mutually_exclusive_group(required=True)
g.add_argument("--out-files", nargs=2, metavar=("DENSITY", "PRESSURE"))
g.add_argument("--mat")
p.add_argument("--cp", type=float, default=662.66)
p.add_argument("--split", type=int, default=None, help="index of top dead centre (default: maximum pressure)")
a = p.parse_args()

if a.mat:
    d = load_mat(a.mat)
    rho, pr, t = (np.ravel(d[k]) for k in ("va_density", "va_pressure", "time"))
else:
    D = np.loadtxt(a.out_files[0], skiprows=3)
    P = np.loadtxt(a.out_files[1], skiprows=3)
    rho, pr, t = D[:, 1], P[:, 1], P[:, 2]
i = int(np.argmax(pr)) if a.split is None else a.split
v = 1 / rho
w_c = trapezoid(pr[: i + 1], v[: i + 1])
w_e = trapezoid(pr[i:], v[i:])
R = 8314 / 28.966
cv = a.cp - R
gamma = a.cp / cv
const = pr[0] * v[0] ** gamma
w_isentropic = quad(lambda x: const / x**gamma, v[0], v[i])[0]
print(f"compression work {w_c:.4g} J/kg, expansion work {w_e:.4g} J/kg")
print(f"isentropic compression work (gamma = {gamma:.4f}) {w_isentropic:.4g} J/kg; closed form "
      f"{const / (1 - gamma) * (v[i] ** (1 - gamma) - v[0] ** (1 - gamma)):.4g} J/kg")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.plot(t, pr, "b", label="expansion")
ax1.plot(t[: i + 1], pr[: i + 1], "r", label="compression")
ax1.set(xlabel="time (s)", ylabel="pressure (Pa)")
ax1.legend()
ax2.plot(v, pr)
ax2.set(xlabel="specific volume (m^3/kg)", ylabel="pressure (Pa)")
plt.show()
