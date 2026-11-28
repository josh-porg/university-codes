"""Gluhareff pressure-jet CFD: build snapshot matrices from the solver's HDF5 output and run DMD
(``data_extracter``, ``perform_DMD_script``, ``DMD_driver`` and the Gluhareff ``DMD`` with
Gavish-Donoho truncation).

    python gluhareff_dmd.py extract Gluhareff_26_2_2025_rawData -o snapshot_matrices
    python gluhareff_dmd.py dmd snapshot_matrices --dt 1e-7 -o DMD_28_2_2025

``extract`` interpolates every field of every ``post*.h5`` file onto the
cell centres of the last file (``scatteredInterpolant`` -> scipy
``LinearNDInterpolator``) and saves ``<FIELD>_snapshot_matrix.npy``.
``dmd`` runs DMD on each matrix and plots the growth-adjusted mode power
against frequency. Needs ``h5py`` (``pip install unicodes[mat73]``).
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import DMD

FIELDS = ["DENSITY", "PRESSURE", "MASSFRAC_O2", "MASSFRAC_H2O", "MASSFRAC_CO", "MASSFRAC_CO2", "MASSFRAC_C3H8",
          "VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z", "VORTICITY_X", "VORTICITY_Y", "VORTICITY_Z", "TEMPERATURE"]

p = argparse.ArgumentParser()
sub = p.add_subparsers(dest="cmd", required=True)
s = sub.add_parser("extract")
s.add_argument("directory", type=Path)
s.add_argument("-o", "--output", type=Path, default=Path("snapshot_matrices"))
s.add_argument("--fields", nargs="+", default=FIELDS)
s = sub.add_parser("dmd")
s.add_argument("directory", type=Path)
s.add_argument("--dt", type=float, default=1e-7)
s.add_argument("--fields", nargs="+", default=["DENSITY", "PRESSURE", "MASSFRAC_C3H8", "VELOCITY_X", "VELOCITY_Y",
                                                "VELOCITY_Z", "VORTICITY_X", "VORTICITY_Y", "VORTICITY_Z",
                                                "TEMPERATURE"])
s.add_argument("-o", "--output", type=Path, default=Path("DMD_results"))
a = p.parse_args()


def cell_centres(f):
    """Mean of each cell's polygon vertices (``cell_centres``)."""
    x, y, z = (np.asarray(f[f"/STREAM_00/VERTEX_COORDINATES/{c}"], float) for c in "XYZ")
    offset = np.asarray(f["/STREAM_00/CONNECTIVITY/POLYGON_OFFSET"], int)
    p2v = np.asarray(f["/STREAM_00/CONNECTIVITY/POLYGON_TO_VERTEX"], int)
    n = len(f["/STREAM_00/CELL_CENTER_DATA/DENSITY"])
    centres = np.empty((n, 3))
    for i in range(n):
        v = p2v[offset[i]:offset[i + 1]]
        centres[i] = x[v].mean(), y[v].mean(), z[v].mean()
    return centres


if a.cmd == "extract":
    import h5py
    from scipy.interpolate import LinearNDInterpolator

    files = sorted(a.directory.glob("post*.h5"), key=lambda q: int(re.search(r"post(\d+)", q.name).group(1)))
    with h5py.File(files[-1]) as f:
        ref = cell_centres(f)
    a.output.mkdir(parents=True, exist_ok=True)
    geometry = []
    for path in files:
        with h5py.File(path) as f:
            geometry.append((cell_centres(f), {k: np.asarray(f[f"/STREAM_00/CELL_CENTER_DATA/{k}"], float)
                                               for k in a.fields}))
    for field in a.fields:
        M = np.column_stack([LinearNDInterpolator(c, d[field])(ref) for c, d in geometry])
        np.save(a.output / f"{field}_snapshot_matrix.npy", M)
        print(f"saved {field}: {M.shape}")
else:
    a.output.mkdir(parents=True, exist_ok=True)
    for field in a.fields:
        X = np.nan_to_num(np.load(a.directory / f"{field}_snapshot_matrix.npy"))
        dmd = DMD(X, dt=a.dt, truncation="svht")
        power = dmd.growth_adjusted_power()
        np.savez(a.output / f"{field}_results.npz", modes=dmd.modes, eigenvalues=dmd.eigenvalues, omega=dmd.omega,
                 time_dynamics=dmd.time_dynamics(), power=power, reconstruction=dmd.reconstruct().real)
        f_hz = np.abs(dmd.frequencies)
        fig, ax = plt.subplots()
        ax.vlines(f_hz, 0, power, colors="b")
        ax.plot(f_hz, power, "ro")
        ax.set(xlabel="Frequency (Hz)", ylabel="Power", title=f"{field} DMD power spectrum (rank {dmd.rank})")
        ax.grid(True)
    plt.show()
