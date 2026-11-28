"""Spectra of signals by FFT and by time-delay (Hankel) DMD (``dmd_psd_playground``, ``PSD_2025``,
``psd_development``)."""

from __future__ import annotations

import numpy as np


def fft_amplitude(x, dt):
    """One-sided FFT amplitude spectrum ``(freqs, |X| 2/N)``."""
    x = np.asarray(x, dtype=float)
    N = len(x)
    X = np.fft.rfft(x)
    return np.fft.rfftfreq(N, dt), np.abs(X) * 2 / N


def hankel(x, rows):
    """Time-delay embedding: row ``k`` is ``x[k : N - rows + k + 1]``."""
    x = np.asarray(x)
    cols = len(x) - rows + 1
    return np.lib.stride_tricks.sliding_window_view(x, cols)[:rows].copy()


def dmd_spectrum(x, dt, delays=500, rank=50):
    """Frequencies (Hz) and amplitudes of the DMD modes of a delay-embedded scalar signal.

    Follows ``dmd_psd_playground``: amplitudes ``|b| 2 / sqrt(delays)`` with
    ``b`` fitted to the first column.
    """
    H = hankel(x, delays)
    U, S, Vh = np.linalg.svd(H[:, :-1], full_matrices=False)
    r = min(rank, len(S))
    U, S, V = U[:, :r], S[:r], Vh[:r].conj().T
    A_t = U.conj().T @ H[:, 1:] @ V / S
    lam, W = np.linalg.eig(A_t)
    Phi = H[:, 1:] @ V / S @ W
    b = np.linalg.lstsq(Phi, H[:, 0], rcond=None)[0]
    freqs = np.log(lam.astype(complex)).imag / (2 * np.pi * dt)
    return freqs, np.abs(b) * 2 / np.sqrt(delays)


def mode_power(time_dynamics):
    """RMS of each mode's time dynamics (``psd_development``), one value per mode (row)."""
    return np.sqrt(np.mean(np.abs(np.asarray(time_dynamics)) ** 2, axis=1))
