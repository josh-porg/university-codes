"""Plot spectra of DMD results exported from PyDMD (``pyDMD Data Analysis``: ``read_data``,
``getModeData``, ``analyze_spectrum``, ``plot_spectrum``, ``data_visualization*``, ``read_modes*``,
``messy_graphs_AAAAAAA``).

The PyDMD run wrote, for each mode ``k``, text files
``<case>_mode_<k>ContTimeEigVal``, ``...Mode`` and ``...TimeDynamic`` holding
Python complex literals, plus ``<case>_initialAmplitudes`` and
``<case>_initialState`` (comma separated)::

    python pydmd_results.py DMD_combustorDataPressure --modes-out modes.npy

Plots: mode power ``|b / omega|`` against frequency ``|Im omega| / 2 pi``
(linear and log axes, with and without the growth/decay scaling),
frequency against mode number and the continuous-time eigenvalues.
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

COMPLEX = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?(?:[-+](?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?j|j)?")


def read_complex(path):
    """All (possibly complex) numbers in a text file written by Python."""
    return np.array([complex(tok) for tok in COMPLEX.findall(Path(path).read_text().replace(" ", ""))])


p = argparse.ArgumentParser()
p.add_argument("case", help="file prefix, e.g. DMD_combustorDataPressure")
p.add_argument("--max-modes", type=int, default=None)
p.add_argument("--modes-out", type=Path, default=None, help="also collect the mode shapes into this .npy")
a = p.parse_args()

omega, modes = [], []
k = 0
while a.max_modes is None or k < a.max_modes:
    f = Path(f"{a.case}_mode_{k}ContTimeEigVal")
    if not f.exists():
        break
    omega.append(read_complex(f)[0])
    if a.modes_out:
        modes.append(read_complex(f"{a.case}_mode_{k}Mode"))
    k += 1
omega = np.array(omega)
b = read_complex(f"{a.case}_initialAmplitudes")[: len(omega)]
print(f"{len(omega)} modes read")
if a.modes_out:
    np.save(a.modes_out, np.column_stack(modes))

freq = np.abs(omega.imag) / (2 * np.pi)
with np.errstate(divide="ignore"):
    power = np.abs(b / omega)
i = np.argmax(power)
print(f"strongest mode {i}: frequency {freq[i]:.4g}, power {power[i]:.4g}")
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
axes[0, 0].plot(freq, power, "bo", ms=3)
axes[0, 0].set(title="power vs frequency (growth/decay scaled)", xlabel="frequency", ylabel="|b / omega|")
axes[0, 1].semilogx(np.sort(freq + 1), power[np.argsort(freq)], "r-x", ms=3)
axes[0, 1].set(title="power vs frequency + 1 (log)", xlabel="frequency + 1")
axes[1, 0].plot(freq, np.abs(b), "bo", ms=3)
axes[1, 0].set(title="initial amplitudes |b| vs frequency", xlabel="frequency")
axes[1, 1].plot(freq, "bo", ms=3)
axes[1, 1].set(title="frequency vs mode number", xlabel="mode number")
fig.tight_layout()
fig, ax = plt.subplots()
ax.scatter(omega.real, omega.imag, s=8)
ax.axhline(0, color="k", lw=0.5)
ax.axvline(0, color="k", lw=0.5)
ax.set(title="continuous-time eigenvalues", xlabel="Re", ylabel="Im")
plt.show()
