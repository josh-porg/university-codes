"""Whole-body vibration ride-quality metrics per ISO 2631-1 (AE 722 ``ride_quality_analysis`` and helpers).

The MATLAB ``design_iso2631_filter`` fed continuous-time polynomial
coefficients straight into the discrete ``filter``; here the ISO 2631-1
Annex A weightings (band limiting, acceleration-velocity transition and
upward step) are built in the s-domain and bilinear-transformed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal
from scipy.integrate import trapezoid

# (f3, f4, Q4, f5, Q5, f6, Q6) per ISO 2631-1 Table A.2; None means no upward step.
_WEIGHTINGS = {
    "Wk": (12.5, 12.5, 0.63, 2.37, 0.91, 3.35, 0.91),
    "Wd": (2.0, 2.0, 0.63, None, None, None, None),
}


def iso2631_weighting(kind: str):
    """Analog weighting filter ``(num, den)`` in s (``design_iso2631_filter``); ``kind`` is 'Wk' (vertical) or 'Wd' (horizontal)."""
    f3, f4, Q4, f5, Q5, f6, Q6 = _WEIGHTINGS[kind]
    w1, w2, w3, w4 = (2 * np.pi * f for f in (0.4, 100.0, f3, f4))
    Q = 1 / np.sqrt(2)
    num = np.polymul([1, 0, 0], [w2**2])  # high-pass x low-pass
    den = np.polymul([1, w1 / Q, w1**2], [1, w2 / Q, w2**2])
    num = np.polymul(num, [1 / w3, 1])  # a-v transition
    den = np.polymul(den, [1 / w4**2, 1 / (Q4 * w4), 1])
    if f5 is not None:
        w5, w6 = 2 * np.pi * f5, 2 * np.pi * f6
        num = np.polymul(num, np.array([1 / w5**2, 1 / (Q5 * w5), 1]) * (w5 / w6) ** 2)
        den = np.polymul(den, [1 / w6**2, 1 / (Q6 * w6), 1])
    return num, den


def apply_weighting(a, fs, kind: str):
    """Frequency-weight an acceleration record sampled at ``fs`` Hz."""
    b, d = signal.bilinear(*iso2631_weighting(kind), fs)
    return signal.lfilter(b, d, a)


def rms(a):
    return float(np.sqrt(np.mean(np.square(a))))


def vibration_dose_value(t, a):
    """VDV = (int a^4 dt)^(1/4) (``calculate_vdv``)."""
    return float(trapezoid(np.asarray(a) ** 4, t) ** 0.25)


def crest_factor(a):
    r = rms(a)
    return float(np.max(np.abs(a)) / r) if r > 0 else np.nan


_RMS_LEVELS = [(0.315, "Not uncomfortable"), (0.63, "A little uncomfortable"), (1.0, "Fairly uncomfortable"),
               (1.6, "Uncomfortable"), (2.5, "Very uncomfortable")]  # fmt: skip
_VDV_LEVELS = [(0.5, "Not uncomfortable"), (1.0, "A little uncomfortable"), (1.4, "Fairly uncomfortable"),
               (2.0, "Uncomfortable"), (2.8, "Very uncomfortable")]  # fmt: skip


def _level(value, table):
    return next((name for limit, name in table if value < limit), "Extremely uncomfortable")


def comfort_level(rms_total, vdv_total):
    """ISO 2631 comfort description; the more severe of the RMS and VDV readings (``assess_comfort``)."""
    by_rms, by_vdv = _level(rms_total, _RMS_LEVELS), _level(vdv_total, _VDV_LEVELS)
    order = [name for _, name in _RMS_LEVELS] + ["Extremely uncomfortable"]
    return max(by_rms, by_vdv, key=order.index)


@dataclass
class RideQuality:
    rms: np.ndarray  # x, y, z
    rms_total: float
    vdv: np.ndarray
    vdv_total: float
    crest: np.ndarray
    comfort: str
    weighted: np.ndarray  # (3, n) weighted accelerations after the transient

    @property
    def figure_of_merit(self) -> float:
        """``rms_total + vdv_total + sum(crest)``, the score compared between rigid and PAH wings."""
        return self.rms_total + self.vdv_total + float(np.sum(self.crest))


def ride_quality(t, ax, ay, az, fs, transient=1.0) -> RideQuality:
    """ISO 2631 analysis of body accelerations (m/s^2) (``ride_quality_analysis``)."""
    w = np.vstack([apply_weighting(ax, fs, "Wd"), apply_weighting(ay, fs, "Wd"), apply_weighting(az, fs, "Wk")])
    k = int(round(transient * fs))
    t, w = np.asarray(t)[k:], w[:, k:]
    r = np.array([rms(c) for c in w])
    v = np.array([vibration_dose_value(t, c) for c in w])
    rt, vt = float(np.sqrt(np.sum(r**2))), float(np.sum(v**4) ** 0.25)
    return RideQuality(r, rt, v, vt, np.array([crest_factor(c) for c in w]), comfort_level(rt, vt), w)
