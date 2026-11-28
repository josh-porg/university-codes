"""Flat-plate boundary layers and skin-friction drag (AE 445, AE 546)."""

from __future__ import annotations

import numpy as np


def cf_laminar(Re):
    """Blasius average skin-friction coefficient ``1.328 / sqrt(Re)`` (``getCfLaminar``)."""
    return 1.328 / np.sqrt(np.asarray(Re, dtype=float))


def cf_turbulent(Re):
    """Prandtl-Schlichting turbulent average skin friction ``0.455 / log10(Re)^2.58``.

    The MATLAB ``getCfTurbulentIncompressable`` used 0.445; 0.455 is the
    textbook constant.
    """
    return 0.455 / np.log10(np.asarray(Re, dtype=float)) ** 2.58


def transition_distance(Re_crit, V, nu):
    """Distance from the leading edge where Re_x reaches ``Re_crit`` (``getLaminarTurbulentBoundryLayerDistance``)."""
    return Re_crit * nu / V


def laminar_thickness(x, Re_x):
    """Blasius boundary-layer thickness ``5 x / sqrt(Re_x)``."""
    return 5 * np.asarray(x) / np.sqrt(Re_x)


def turbulent_thickness(x, Re_x):
    """1/7-power-law turbulent thickness ``0.37 x / Re_x^(1/5)``."""
    return 0.37 * np.asarray(x) / np.asarray(Re_x) ** 0.2


def flat_plate_drag(V, chord, span, rho, nu, Re_crit=None, sides=2):
    """Skin-friction drag of a flat plate at zero incidence (AE 445 HW 7).

    Fully turbulent if ``Re_crit`` is None, otherwise laminar up to the
    transition point and turbulent after it (turbulent part computed as the
    full-plate turbulent drag minus a turbulent run of the laminar length).
    """
    q = 0.5 * rho * V**2
    Re_L = V * chord / nu
    if Re_crit is None:
        return sides * q * cf_turbulent(Re_L) * chord * span
    x_t = min(transition_distance(Re_crit, V, nu), chord)
    laminar = q * cf_laminar(Re_crit if x_t < chord else Re_L) * x_t * span
    if x_t >= chord:
        return sides * laminar
    turbulent = q * span * (cf_turbulent(Re_L) * chord - cf_turbulent(Re_crit) * x_t)
    return sides * (laminar + turbulent)
