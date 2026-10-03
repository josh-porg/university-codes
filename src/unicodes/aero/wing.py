"""Wing aerodynamics: lift-curve slopes, induced drag, sweep and compressibility (AE 445, AE 550, AE 722).

Angles in radians. The MATLAB names are given in each docstring.
"""

from __future__ import annotations

import numpy as np


# Lift-curve slope


def finite_wing_lift_slope(a0, A, e=1.0):
    """Lift-curve slope of a finite wing from its airfoil slope ``a0`` (``getCLAlphaWIngFromAirfoil``).

    ``a = a0 / (1 + a0 / (pi A e))``; good for A > 5 and little sweep.
    """
    return a0 / (1 + a0 / (np.pi * A * e))


def rectangular_wing_lift_slope(A):
    """Incompressible lift-curve slope of a rectangular wing (``getCLaplhaRectangularIncompresable``)."""
    return 2 * np.pi * A / (2 + np.sqrt(A**2 + 4))


def helmbold_lift_slope(A, M=0.0, sweep_half_chord=0.0, a0=2 * np.pi):
    """DATCOM/Helmbold lift-curve slope for a swept wing in subsonic flow (``getCLalpha``).

    ``CL_a = 2 pi A / (2 + sqrt(A^2 beta^2 / k^2 (1 + tan^2 L_c/2 / beta^2) + 4))``
    with ``beta = sqrt(1 - M^2)`` and ``k = a0 / (2 pi)``. (The MATLAB version
    used ``beta = 1 - M^2``, divided ``k`` by beta and referenced an undefined ``B``.)
    """
    beta2 = 1 - np.asarray(M, dtype=float) ** 2
    k = a0 / (2 * np.pi)
    return 2 * np.pi * A / (2 + np.sqrt(A**2 * beta2 / k**2 * (1 + np.tan(sweep_half_chord) ** 2 / beta2) + 4))


def rescale_lift_slope(a1, A_eff1, A_eff2):
    """Lift slope of a second wing with the same airfoil, e.g. model to full scale (``getCLalphafromClAlpha``)."""
    return a1 / (1 + a1 / np.pi * (1 / A_eff2 - 1 / A_eff1))


def sweep_at(x_c, sweep_le, A, taper):
    """Sweep of the ``x_c`` chord line (e.g. 0.25, 0.5) of a straight-tapered wing from the LE sweep."""
    return np.arctan(np.tan(sweep_le) - 4 * x_c * (1 - taper) / (A * (1 + taper)))


def polhamus_lift_slope(a0, A, sweep_le, taper, M=0.0):
    """Polhamus lift-curve slope of a planform (``polhamus``).

    Better than Prandtl-Glauert in the transonic range but still not valid at M = 1.
    ``k`` uses the original course correlation with ``sweep_le`` in radians.
    """
    tan_half = np.tan(sweep_at(0.5, sweep_le, A, taper))
    if A < 4:
        k = 1 + A * (1.87 - 0.000233 * sweep_le) / 100
    else:
        k = 1 + ((8.2 - 2.3 * sweep_le) - A * (0.22 - 0.153 * sweep_le)) / 100
    beta2 = 1 - np.asarray(M, dtype=float) ** 2
    return a0 * A / (2 + np.sqrt(A**2 * beta2 / k**2 * (1 + tan_half**2 / beta2) + 4))


# Drag polar


def effective_aspect_ratio(A, e):
    """``A e`` (``getEffectiveAspectRatio``)."""
    return A * e


def induced_drag_coefficient(C_L, A, e=1.0):
    """``C_Di = C_L^2 / (pi A e)`` (``getInducedDragCoefficient``)."""
    return np.asarray(C_L) ** 2 / (np.pi * A * e)


def induced_angle(C_L, A, e=1.0):
    """Induced angle of attack ``C_L / (pi A e)`` (rad).

    The MATLAB ``getInducedAngleOfAttack`` returned ``C_L / C_Di``, which is not an angle.
    """
    return np.asarray(C_L) / (np.pi * A * e)


def parabolic_drag_polar(C_L, C_D0, A, e):
    """``C_D = C_D0 + C_L^2 / (pi A e)`` (``DragPolar_parabolic``, ``getParabolicDragPolar``)."""
    return C_D0 + induced_drag_coefficient(C_L, A, e)


def max_lift_to_drag(C_D0, A, e):
    """``(L/D)_max = 0.5 sqrt(pi A e / C_D0)``."""
    return 0.5 * np.sqrt(np.pi * A * e / C_D0)


def lift_to_drag(wing_loading, rho, V, C_D0, A, e):
    """L/D in level flight at speed ``V`` (``LD``)."""
    C_L = 2 * wing_loading / (rho * np.asarray(V, dtype=float) ** 2)
    return C_L / parabolic_drag_polar(C_L, C_D0, A, e)


def stall_speed(wing_loading, C_L_max, rho):
    """1 g stall speed ``sqrt(2 W/S / (rho C_Lmax))`` (``StallSpeed``, ``getStallSpeed``)."""
    return np.sqrt(2 * np.asarray(wing_loading) / (rho * C_L_max))


# Compressibility


def prandtl_glauert(value_incompressible, M):
    """Prandtl-Glauert correction ``value / sqrt(1 - M^2)``."""
    return value_incompressible / np.sqrt(1 - np.asarray(M, dtype=float) ** 2)


def karman_tsien(Cp0, M):
    """Karman-Tsien compressibility correction of a pressure coefficient."""
    M = np.asarray(M, dtype=float)
    b = np.sqrt(1 - M**2)
    return Cp0 / (b + Cp0 * M**2 / (2 * (1 + b)))


def laitone(Cp0, M, gamma=1.4):
    """Laitone compressibility correction of a pressure coefficient."""
    M = np.asarray(M, dtype=float)
    b = np.sqrt(1 - M**2)
    return Cp0 / (b + Cp0 * M**2 * (1 + (gamma - 1) / 2 * M**2) / (2 * (1 + b)))


def critical_mach_swept_2d(M_crit_unswept, sweep_le):
    """Critical Mach of an infinite swept wing ``M* / cos(L)`` (``getCriticalMachNumber2D``)."""
    return M_crit_unswept / np.cos(sweep_le)


def critical_mach_swept_3d(M_crit_unswept, sweep_le):
    """Critical Mach of a finite swept wing ``M* / sqrt(cos L)`` (``getCriticalMachNumber3D``)."""
    return M_crit_unswept / np.sqrt(np.cos(sweep_le))


def effective_mach(M, sweep_le):
    """Mach number normal to a swept leading edge (``getMachNumberEffective``)."""
    return M * np.cos(sweep_le)


def critical_mach_from_cp(Cp0, method: str = "karman_tsien", gamma=1.4):
    """Free-stream Mach at which the corrected minimum Cp reaches the sonic value (AE 445 HW 12)."""
    from scipy.optimize import brentq

    from ..gasdynamics import critical_pressure_coefficient

    correct = {"prandtl_glauert": prandtl_glauert, "karman_tsien": karman_tsien}[method]
    return brentq(lambda M: correct(Cp0, M) - critical_pressure_coefficient(M, gamma), 1e-3, 0.999)


# Planform geometry


def mean_geometric_chord(c_r, b, taper, sweep_le=0.0):
    """Mean geometric chord of a straight-tapered panel and its position (``MeanGeometricCord``).

    Returns ``(c_bar, x_mgc, y_mgc)``: length, chordwise offset of its leading
    edge from the root LE, and spanwise station.
    """
    c_bar = 2 / 3 * c_r * (1 + taper + taper**2) / (1 + taper)
    y = b * (1 + 2 * taper) / (6 * (1 + taper))
    return c_bar, y * np.tan(sweep_le), y


def root_chord_from_mgc(c_bar, taper):
    """Root chord of a straight-tapered wing with mean geometric chord ``c_bar``."""
    return 1.5 * c_bar * (1 + taper) / (taper**2 + taper + 1)


def schrenk_lift_distribution(y, b, c_r, c_t, total_lift):
    """Schrenk's lift per unit span: average of elliptic and trapezoidal loadings (``Schrenk_Lift_aproximation``).

    ``y`` is the spanwise station from the root (0 to b/2); ``total_lift``
    is for the whole wing.
    """
    y = np.asarray(y, dtype=float)
    eta = 2 * y / b
    elliptic = 4 * total_lift / (np.pi * b) * np.sqrt(np.clip(1 - eta**2, 0, None))
    trapezoidal = 2 * total_lift / (b * (1 + c_t / c_r)) * (1 - eta * (1 - c_t / c_r))
    return (elliptic + trapezoidal) / 2


# Longitudinal stability (AE 550, AE 722)


def downwash_gradient(m, r, A, sweep_le, taper, C_L_alpha_ratio=1.0):
    """DATCOM downwash gradient d(eps)/d(alpha) at the tail (``downwashGradient``).

    ``r``, ``m``: horizontal and vertical distances between the wing and tail
    aerodynamic centres divided by the wing half span. ``C_L_alpha_ratio`` is
    ``C_L_alpha(M) / C_L_alpha(0)`` for the compressibility correction.
    """
    sweep_qc = sweep_at(0.25, sweep_le, A, taper)
    K_A = 1 / A - 1 / (1 + A**1.7)
    K_lambda = (10 - 3 * taper) / 7
    K_mr = (1 - m / 2) / r**0.33
    return 4.44 * (K_A * K_lambda * K_mr * np.sqrt(np.cos(sweep_qc))) ** 1.19 * C_L_alpha_ratio


def cm_alpha_wing_fuselage(C_L_alpha_wf, x_cg, x_ac_wf):
    """Wing-fuselage pitch stiffness ``(x_cg - x_ac) C_L_alpha_wf`` (positions over c_bar)."""
    return (x_cg - x_ac_wf) * C_L_alpha_wf


def cm_alpha_tail(S_h, S_ref, eta_h, C_L_alpha_h, downwash_grad, x_cg, x_ac_h):
    """Horizontal-tail contribution to C_m_alpha."""
    return S_h / S_ref * eta_h * C_L_alpha_h * (1 - downwash_grad) * (x_cg - x_ac_h)


def cm_tail_incidence(S_h, S_ref, eta_h, C_L_alpha_h, x_cg, x_ac_h):
    """C_m due to horizontal-tail incidence (positive trailing edge up)."""
    return -S_h / S_ref * eta_h * C_L_alpha_h * (x_cg - x_ac_h)


def cm0_tail(S_h, S_ref, eta_h, C_L_alpha_h, epsilon_0, x_cg, x_ac_h):
    """Horizontal-tail C_m at zero angle of attack from the zero-alpha downwash."""
    return C_L_alpha_h * eta_h * S_h / S_ref * epsilon_0 * (x_cg - x_ac_h)


# Free vortex (AE 445 HW 4)


def free_vortex_velocity(r, circulation):
    """Tangential velocity of a free vortex ``-Gamma / (2 pi r)``."""
    return -circulation / (2 * np.pi * np.asarray(r, dtype=float))
