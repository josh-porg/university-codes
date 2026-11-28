"""POD and DMD of the turbulent separation bubble PIV data
(``Proper_orthogonal_decomposition_turbulent_seperation_bubble``, ``DMDTurbulentSeperationBubble``).

Needs ``data_tsb.mat`` (``U_Z0_VER_Run1``, ``V_Z0_VER_Run1``, ``axeX_Z0``,
``axeY_Z0``); the fields are subsampled by 3 in space as in the MATLAB::

    python separation_bubble.py data_tsb.mat --method pod --modes 5
    python separation_bubble.py data_tsb.mat --method dmd --dt 1.1 --show-modes 6

``--animate out.gif`` writes the reconstruction next to the data.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from unicodes.decomposition import DMD, POD
from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("data")
p.add_argument("--method", choices=["pod", "dmd"], default="pod")
p.add_argument("--component", choices=["U", "V"], default="U")
p.add_argument("--modes", type=int, default=5, help="modes in the reconstruction")
p.add_argument("--show-modes", type=int, default=4)
p.add_argument("--dt", type=float, default=1.1)
p.add_argument("--frames", type=int, default=56)
p.add_argument("--animate", default=None)
a = p.parse_args()

d = load_mat(a.data)
X = np.asarray(d["axeX_Z0"]).T[::3, ::3]
Yr = np.asarray(d["axeY_Z0"])
Y = np.flipud(Yr - Yr[0, 0]).T[::3, ::3]
field = np.transpose(np.asarray(d[f"{a.component}_Z0_VER_Run1"], float), (1, 0, 2))[::3, ::3, :]
field = np.nan_to_num(field)

if a.method == "pod":
    model = POD(field)
    recon = model.reconstruct(a.modes)
    modes = model.get_modes()[..., : a.show_modes]
    print("energy fraction of the first modes:", np.round(model.energy_fraction[:10], 4))
else:
    model = DMD(field, dt=a.dt, truncation=0.005)
    idx = np.argsort(-np.abs(model.amplitudes))[: a.modes]
    recon = model._to_data(model.reconstruct(modes=idx).real)
    modes = model.get_modes()[..., np.argsort(-np.abs(model.amplitudes))[: a.show_modes]].real
    print(f"DMD rank {model.rank}; leading frequencies:", np.round(np.abs(model.frequencies[idx]), 4))
    fig, ax = plt.subplots()
    ax.plot(model.omega.real, model.omega.imag, "+r")
    ax.set(title="continuous-time eigenvalues")

fig, axes = plt.subplots(int(np.ceil(modes.shape[-1] / 2)), 2, figsize=(11, 2.5 * modes.shape[-1] / 2 + 1))
for k, ax in enumerate(np.ravel(axes)[: modes.shape[-1]]):
    ax.contourf(X, Y, modes[..., k], 30, cmap="RdBu_r")
    ax.set_title(f"{a.method.upper()} mode {k + 1}")
fig.tight_layout()

fig, (ax1, ax2) = plt.subplots(2, 1)
levels = np.linspace(-7, 20, 30)


def frame(n):
    for ax, data, title in ((ax1, recon, f"{a.modes}-mode reconstruction"), (ax2, field, "data")):
        ax.clear()
        ax.contourf(X, Y, data[..., n], levels, extend="both")
        ax.set_title(f"{title}, snapshot {n}")


frame(0)
if a.animate:
    FuncAnimation(fig, frame, frames=min(a.frames, field.shape[-1]), interval=50).save(a.animate)
plt.show()
