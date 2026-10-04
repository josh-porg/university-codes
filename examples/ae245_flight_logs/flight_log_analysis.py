"""AE 245: ArduPilot flight-log analysis (``ae_245_analysis``, ``ae245quickanalysis``, ``ae_245_hw_3``,
``lab2_1``).

Pass one or more ``.BIN-*.mat`` logs converted by Mission Planner (not in the
repository), e.g. ``"Maiden Voyage.BIN-29424.mat"``::

    python flight_log_analysis.py "Maiden Voyage.BIN-29424.mat" "Roll Autotune (12-9).BIN-29538.mat"

For each log: attitude tracking errors (ATT demanded vs achieved) with the
descriptive statistics of the MATLAB (RMS, mean, median, std, MAD, IQR,
skewness, kurtosis, correlation, PCA explained variance), altitude and
climb-rate tracking (CTUN), body-rate tracking (RATE), the flat-earth GPS
path (``lla2flat``; the MATLAB had this switched off as broken) and the
HW 3 position/attitude panels.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from unicodes.flight_dynamics.kinematics import lla_to_flat
from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("logs", nargs="+")
a = p.parse_args()


def describe(name, demanded, actual):
    e = demanded - actual
    q75, q25 = np.percentile(e, [75, 25])
    pca = np.linalg.svd(np.column_stack([actual, demanded]) - np.mean([actual, demanded], axis=1), compute_uv=False) ** 2
    print(f"  {name:6s} RMS {np.sqrt(np.mean(e**2)):8.3f}  mean {e.mean():8.3f}  median {np.median(e):8.3f}  "
          f"std {e.std(ddof=1):8.3f}  MAD {np.mean(np.abs(e - e.mean())):7.3f}  IQR {q75 - q25:7.3f}  "
          f"skew {stats.skew(e):6.2f}  kurt {stats.kurtosis(e, fisher=False):6.2f}  "
          f"r {np.corrcoef(actual, demanded)[0, 1]:6.3f}  PCA {100 * pca[0] / pca.sum():5.1f} %")


for log in a.logs:
    d = load_mat(log)
    name = Path(log).name.split(".")[0]
    print(name)
    ATT = np.asarray(d["ATT"])
    t = ATT[:, 1]
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 7))
    for ax, (label, i_d, i_a) in zip(axes, (("roll", 2, 3), ("pitch", 4, 5), ("yaw", 6, 7))):
        describe(label, ATT[:, i_d], ATT[:, i_a])
        ax.plot(t, ATT[:, i_d], label="demanded")
        ax.plot(t, ATT[:, i_a], label="achieved")
        ax.set_ylabel(f"{label} (deg)")
    axes[0].legend()
    axes[0].set_title(f"{name}: attitude tracking")
    if "CTUN" in d:
        C = np.asarray(d["CTUN"])
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.plot(C[:, 1], C[:, 6], label="demanded")
        ax1.plot(C[:, 1], C[:, 7], label="achieved")
        ax1.set_ylabel("altitude")
        ax2.plot(C[:, 1], C[:, 12], C[:, 1], C[:, 13])
        ax2.set_ylabel("climb rate")
        ax1.legend()
        ax1.set_title(f"{name}: altitude control")
    if "RATE" in d:
        Rt = np.asarray(d["RATE"])
        fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8, 8))
        for ax, (label, i_d, i_a) in zip(axes, (("p", 2, 3), ("q", 5, 6), ("r", 8, 9), ("a_z", 11, 12))):
            ax.plot(Rt[:, 1], Rt[:, i_d], Rt[:, 1], Rt[:, i_a])
            ax.set_ylabel(label)
        axes[0].set_title(f"{name}: rate tracking (desired, achieved)")
    if "GPS" in d:
        G = np.asarray(d["GPS"])
        n, e, dn = lla_to_flat(G[:, 7], G[:, 8], G[:, 9], G[0, 7], G[0, 8], -G[0, 9])
        fig = plt.figure(figsize=(11, 7))
        ax = fig.add_subplot(1, 2, 1, projection="3d")
        ax.plot(e, n, -dn)
        ax.set(xlabel="east (m)", ylabel="north (m)", zlabel="up (m)", title=f"{name}: GPS path")
        for k, (lab, col) in enumerate((("roll", 3), ("pitch", 5), ("yaw", 7))):
            axk = fig.add_subplot(3, 2, 2 * k + 2)
            axk.plot(t, ATT[:, col])
            axk.set_ylabel(lab)
plt.show()
