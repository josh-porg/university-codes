"""AE 430 instrumentation: moving averages, harmonic identities, FFT spectra and filtered-signal
reconstruction (``HW 1/Problem_1``, ``HW 2/HW_2``, ``Lab 7/Lab_7_Data_analysis_version_1``).

Data files are not in the repository; each part runs when its file is given::

    python signal_analysis.py --noisy noisy.txt --signal signal.mat --lab7-dir "Lab 7"

Lab 7 reads the LabVIEW ``.lvm`` files (two columns: raw, filtered) for
60/80/100 Hz square and triangle waves through low- and high-pass filters,
sampled at 2 kHz.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import fft_amplitude
from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("--noisy")
p.add_argument("--signal")
p.add_argument("--lab7-dir", type=Path)
a = p.parse_args()

# HW 1: 5 sin 4t + 3 cos 4t = R sin(4t + phi)
R, phi = np.hypot(5, 3), np.arctan2(3, 5)
print(f"HW 1: 5 sin(4t) + 3 cos(4t) = {R:.4f} sin(4t + {phi:.4f}) = {R:.4f} cos(4t - {np.pi / 2 - phi:.4f})")
if a.noisy:
    x = np.loadtxt(a.noisy).ravel()
    fig, ax = plt.subplots()
    ax.plot(x, label="raw data")
    for n in (2, 3, 4):  # movmean: centred window, shrinking at the ends
        ax.plot([x[max(0, i - n // 2):i + (n - 1) // 2 + 1].mean() for i in range(len(x))], label=f"{n}-point moving average")
    ax.set(title="Amplitude vs time", xlabel="Time", ylabel="Amplitude")
    ax.legend()

if a.signal:
    d = load_mat(a.signal)
    X = np.ravel(d["signalB6"])
    Fs = 1000
    f, P1 = fft_amplitude(X, 1 / Fs)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(Fs / len(X) * (np.arange(len(X)) - len(X) // 2), np.abs(np.fft.fftshift(np.fft.fft(X))), lw=2)
    ax1.set(title="FFT spectrum, positive and negative frequencies", xlabel="f (Hz)", ylabel="|fft(X)|")
    ax2.plot(f, P1)
    ax2.set(title="single-sided amplitude spectrum", xlabel="f (Hz)")
    print("HW 2: dominant frequencies (Hz):", np.round(f[np.argsort(P1)[-3:]], 1))

if a.lab7_dir:
    fs, n = 2000.0, 2000
    t = np.arange(n) / fs
    for path in sorted(a.lab7_dir.glob("*.lvm")):
        data = np.loadtxt(path)[:n, :2]
        f, amp_raw = fft_amplitude(data[:, 0], 1 / fs)
        _, amp_filt = fft_amplitude(data[:, 1], 1 / fs)
        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.plot(t[:50], data[:50, 0], label="raw")
        ax1.plot(t[:50], data[:50, 1], label="filtered")
        ax1.set(title=f"{path.stem} recorded signals", xlabel="t (s)", ylabel="V")
        ax1.legend()
        ax2.plot(f, amp_raw, label="raw")
        ax2.plot(f, amp_filt, label="filtered")
        ax2.set(title="amplitude spectrum", xlabel="f (Hz)")
        ax2.legend()
plt.show()
