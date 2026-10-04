"""DMD with control on flight-test data (``DMDc/DMDc``).

Needs ``FlightDataUncontrolled.mat`` holding ``StateAndControl`` (five state
rows, one control row)::

    python dmdc_flight.py FlightDataUncontrolled.mat

The identified continuous-time eigenvalues are compared with those of the
longitudinal model used in the MATLAB script.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import dmdc
from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("data")
p.add_argument("--dt", type=float, default=1.0)
p.add_argument("--tol", type=float, default=1e-5)
a = p.parse_args()

SC = np.asarray(load_mat(a.data)["StateAndControl"], dtype=float)
res = dmdc(SC[:5], SC[5:6], dt=a.dt, tol=a.tol)
print("DMDc continuous-time eigenvalues:", res.omega)
A_long = np.array([[-0.0374, 17.4632, 0, -32.0236], [-0.340, -170.4133, 185.0700, -3.3658],
                   [0, -5.5002, -2.1809, 0], [0, 0, 1.0, 0]])
ref = np.linalg.eigvals(A_long)
print("Reference longitudinal model eigenvalues:", ref)
fig, ax = plt.subplots()
ax.plot(res.omega.real, res.omega.imag, "xr", label="DMDc")
ax.plot(ref.real, ref.imag, "xk", label="model")
ax.set(xlabel="Re", ylabel="Im")
ax.legend()
plt.show()
