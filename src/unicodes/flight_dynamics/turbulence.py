"""Dryden turbulence per MIL-HDBK-1797 (AE 722 ``generate_dryden_gust``, ``..._with_rotational``).

Shaping filters are discretised with the bilinear (Tustin) transform and
driven by white noise, as in the MATLAB. Velocities in m/s.

The MATLAB drove the filters with unit-variance ``randn``; for the MIL
spectra (written with the ``1/pi`` factor) the discrete noise needs variance
``pi/dt`` for the output RMS to equal sigma, so the old gusts were about
``sqrt(dt/pi)`` too weak (0.18x at 10 Hz). The noise is scaled here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

FT = 0.3048
KT = 0.514444


@dataclass
class DrydenParameters:
    sigma_u: float  # m/s
    sigma_v: float
    sigma_w: float
    L_u: float  # m
    L_v: float
    L_w: float


def dryden_parameters(altitude, intensity="moderate", sigma_high=15 * FT):
    """Turbulence intensities and scale lengths for an altitude (m) above ground.

    Below 1000 ft the low-altitude model with the 20 ft wind for
    ``intensity`` (15/30/45 kt); above 2000 ft isotropic with
    ``L_u = 2 L_v = L_w = 1750 ft`` and ``sigma = sigma_high`` (the MATLAB
    hard-coded the severe 46,000 ft value, 15 ft/s); linear blend between.
    """
    W20 = {"light": 15, "moderate": 30, "severe": 45}[intensity.lower()] * KT
    h = altitude / FT

    def low(h):
        sw = 0.1 * W20
        su = sw / (0.177 + 0.000823 * h) ** 0.4
        Lu = h / (0.177 + 0.000823 * h) ** 1.2
        return np.array([su, su, sw, Lu, Lu / 2, h])

    high = np.array([sigma_high, sigma_high, sigma_high, 1750, 875, 1750])
    if h < 1000:
        p = low(max(h, 10.0))
    elif h > 2000:
        p = high
    else:
        a = (h - 1000) / 1000
        p = low(1000.0) + a * (high - low(1000.0))
    return DrydenParameters(*p[:3], *(p[3:] * FT))


def _filter(num, den, dt, noise):
    b, a, _ = signal.cont2discrete((num, den), dt, method="bilinear")
    return signal.lfilter(np.ravel(b), a, noise * np.sqrt(np.pi / dt))


def dryden_gusts(n, dt, V, altitude, intensity="moderate", rng=None, sigma_high=15 * FT):
    """Linear gust components ``(u_g, v_g, w_g)`` (m/s), ``n`` samples at ``dt`` (MIL spec suggests 10 Hz)."""
    rng = np.random.default_rng(rng)
    p = dryden_parameters(altitude, intensity, sigma_high)
    Tu, Tv, Tw = p.L_u / V, p.L_v / V, p.L_w / V
    ug = _filter([p.sigma_u * np.sqrt(2 * Tu / np.pi)], [Tu, 1], dt, rng.standard_normal(n))

    def lateral(sigma, T):
        k = sigma * np.sqrt(2 * T / np.pi)
        return [k * 2 * np.sqrt(3) * T, k], np.polymul([2 * T, 1], [2 * T, 1])

    vg = _filter(*lateral(p.sigma_v, Tv), dt, rng.standard_normal(n))
    wg = _filter(*lateral(p.sigma_w, Tw), dt, rng.standard_normal(n))
    return ug, vg, wg


def dryden_rotational_gusts(wg, vg, dt, V, altitude, semi_span, rng=None):
    """Rotational gusts ``(p_g, q_g, r_g)`` (rad/s) consistent with ``w_g``/``v_g`` (``generate_dryden_gust_with_rotational``).

    As in the MATLAB, sigma_w is taken from the RMS of ``wg`` and ``q_g``,
    ``r_g`` are driven by fresh noise through ``G_w``/``G_v`` shaped filters.
    """
    rng = np.random.default_rng(rng)
    h = altitude / FT
    Lw = (h if h < 1000 else 1750 if h > 2000 else np.interp(h, [1000, 2000], [1000, 1750])) * FT
    sigma = np.std(wg)
    n = len(wg)
    b = semi_span
    tau_p = 4 * b / (np.pi * V)
    K = sigma * np.sqrt(0.8 / V) * (np.pi / (4 * b)) ** (1 / 6) / (2 * Lw) ** (1 / 3)
    pg = _filter([K], [tau_p, 1], dt, rng.standard_normal(n))

    def Gw(L):
        T = L / V
        k = sigma * np.sqrt(2 * T / np.pi)
        return np.array([k * 2 * np.sqrt(3) * T, k]), np.polymul([2 * T, 1], [2 * T, 1])

    num, den = Gw(Lw)
    qg = _filter(np.polymul([1 / V, 0], num), np.polymul([tau_p, 1], den), dt, rng.standard_normal(n))
    num, den = Gw(Lw / 2)
    tau_r = 3 * b / (np.pi * V)
    rg = _filter(-np.polymul([1 / V, 0], num), np.polymul([tau_r, 1], den), dt, rng.standard_normal(n))
    return pg, qg, rg
