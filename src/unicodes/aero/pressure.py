"""Section force coefficients from surface pressure taps and manometer readings (AE 546 labs)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

G0 = 9.81


def manometer_velocity(h_total, h_static, rho_liquid, rho_air, g=G0):
    """Free-stream speed from a pitot-static manometer head difference (same length units -> m)."""
    return np.sqrt(2 * rho_liquid * g * np.abs(np.asarray(h_total) - h_static) / rho_air)


def head_to_pressure(h, rho_liquid, g=G0):
    """Gauge pressure of a liquid column."""
    return rho_liquid * g * np.asarray(h)


def cp_from_heads(h_local, h_static, h_total):
    """C_p from manometer heads, ``(h - h_s) / (h_s - h_t)``.

    Manometer heads fall as pressure rises, so the sign is flipped relative to
    :func:`unicodes.gasdynamics.pressure_coefficient`.
    """
    return (np.asarray(h_local) - h_static) / (h_static - h_total)


def tap_geometry(x_surface, z_surface):
    """Tap locations (midpoints) and panel extents ``dx, dz`` along a surface polyline."""
    x, z = np.asarray(x_surface, dtype=float), np.asarray(z_surface, dtype=float)
    return (x[:-1] + x[1:]) / 2, (z[:-1] + z[1:]) / 2, np.diff(x), np.diff(z)


@dataclass
class SectionCoefficients:
    c_n: np.ndarray
    c_a: np.ndarray
    c_l: np.ndarray
    c_d: np.ndarray
    c_m_le: np.ndarray
    c_m_qc: np.ndarray


def section_coefficients(alpha, cp_upper, cp_lower, upper, lower) -> SectionCoefficients:
    """Integrate tap pressures into normal/axial, lift/drag and moment coefficients.

    ``upper``/``lower`` are ``(x, z)`` polylines (chord-normalised) running
    leading edge to trailing edge; ``cp_upper``/``cp_lower`` hold one value per
    panel (rows) and may have one column per angle of attack ``alpha`` (rad).
    Moments are positive nose-up.
    """
    xu, zu, dxu, dzu = tap_geometry(*upper)
    xl, zl, dxl, dzl = tap_geometry(*lower)
    cpu, cpl = np.asarray(cp_upper, dtype=float), np.asarray(cp_lower, dtype=float)
    c_n = -(dxu @ cpu) + dxl @ cpl
    c_a = dzu @ cpu - dzl @ cpl
    c_m_le = (xu * dxu) @ cpu - (xl * dxl) @ cpl + (zu * dzu) @ cpu - (zl * dzl) @ cpl
    alpha = np.asarray(alpha, dtype=float)
    c_l = c_n * np.cos(alpha) - c_a * np.sin(alpha)
    c_d = c_n * np.sin(alpha) + c_a * np.cos(alpha)
    return SectionCoefficients(c_n, c_a, c_l, c_d, c_m_le, c_m_le + c_l / 4)
