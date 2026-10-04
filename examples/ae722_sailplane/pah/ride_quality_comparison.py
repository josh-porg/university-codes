"""ISO 2631 ride quality of rigid vs PAH aircraft (``ride_quality_analysis``, ``ride_quality_analsyis(_V1)``,
``load_example_data``, ``plot_ISO2631_results``, ``display_ISO2631_results``).

With no arguments uses synthetic flight data like ``load_example_data``;
otherwise pass ``.mat`` files holding a ``flight_data`` struct
(``t, ax, ay, az, fs``) for the control and the test aircraft.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from unicodes.flight_dynamics.ride_quality import ride_quality


def example_data(fs=100.0, T=60.0, rng=0):
    rng = np.random.default_rng(rng)
    t = np.arange(0, T, 1 / fs)
    az = 0.5 * np.sin(2 * np.pi * 1.5 * t) + 0.3 * rng.standard_normal(t.size)
    ax = 0.1 * rng.standard_normal(t.size)
    ay = 0.15 * rng.standard_normal(t.size)
    return dict(t=t, ax=ax, ay=ay, az=az, fs=fs)


def load(path):
    from unicodes.io import load_mat

    d = load_mat(path)["flight_data"]
    return {k: np.ravel(d[k]) if k != "fs" else float(np.ravel(d[k])[0]) for k in ("t", "ax", "ay", "az", "fs")}


runs = {"control": example_data(rng=0), "test": example_data(rng=1)}
if len(sys.argv) == 3:
    runs = {"control": load(sys.argv[1]), "test": load(sys.argv[2])}

for name, d in runs.items():
    r = ride_quality(d["t"], d["ax"], d["ay"], d["az"], d["fs"])
    print(f"{name:8s} RMS {r.rms_total:.3f} m/s^2, VDV {r.vdv_total:.3f}, crest {np.round(r.crest, 2)}, "
          f"'{r.comfort}', FOM {r.figure_of_merit:.3f}")
    f, p_raw = signal.welch(d["az"], d["fs"], nperseg=1024)
    f, p_w = signal.welch(r.weighted[2], d["fs"], nperseg=1024)
    plt.semilogx(f, 10 * np.log10(p_raw), label=f"{name} raw")
    plt.semilogx(f, 10 * np.log10(p_w), "--", label=f"{name} weighted (Wk)")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Vertical acceleration PSD (dB/Hz)")
plt.legend()
plt.show()
