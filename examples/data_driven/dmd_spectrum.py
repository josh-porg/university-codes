"""Power spectra from FFT and from time-delay DMD (``dmd_psd_playground``, ``PSD_2025``).

A noisy two-tone signal (7 Hz and 13 Hz) and a clean three-tone signal
(5, 50, 120 Hz) are analysed with both methods. ``PSD_2025`` stacked
circularly shifted copies of the signal and took ``|Phi|^2`` of the first
mode row as power; the delay-embedding amplitudes of the playground script
are used for both here.
"""

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import dmd_spectrum, fft_amplitude

rng = np.random.default_rng(0)
cases = []
dt = 0.01
t = np.arange(0, 10 + dt / 2, dt)
clean = 14 * np.sin(7 * 2 * np.pi * t) + 5 * np.sin(13 * 2 * np.pi * t)
cases.append(("Noisy 7 + 13 Hz", t, clean + 10 * rng.standard_normal(t.size), clean, 500, 50))
t2 = np.linspace(0, 1, 1000)
sig = sum(A * np.sin(2 * np.pi * f * t2) for f, A in zip((5, 50, 120), (1, 0.5, 0.2)))
cases.append(("5 + 50 + 120 Hz", t2, sig, sig, 100, 10))

for title, tt, x, xc, delays, rank in cases:
    h = tt[1] - tt[0]
    f_fft, p_fft = fft_amplitude(x, h)
    f_dmd, p_dmd = dmd_spectrum(x, h, delays=delays, rank=rank)
    keep = f_dmd >= 0
    print(title, "DMD peaks (Hz):", np.round(np.sort(f_dmd[keep][np.argsort(p_dmd[keep])[-3:]]), 2))
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(tt, x, color=(0.7, 0.3, 0.3), lw=0.9, label="signal")
    ax1.plot(tt, xc, "k", lw=1.5, label="clean")
    ax1.legend()
    ax2.plot(f_fft, p_fft, "k", lw=1.2, label="FFT")
    ax2.scatter(np.abs(f_dmd), p_dmd, c="r", label="DMD")
    ax2.set(xlabel="Frequency (Hz)", ylabel="Amplitude", title=title)
    ax2.legend()
plt.show()
