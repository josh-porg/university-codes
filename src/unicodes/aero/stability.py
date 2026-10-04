"""Static stability build-ups for wing-body + horizontal tail/canard + vertical tail (AE 550 exams).

Positions are fractions of the wing mean geometric chord (x positive aft);
``ac_shift`` is a Munk/body shift of the wing-body aerodynamic centre.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Surface:
    """A horizontal lifting surface behind (tail) or ahead of (canard) the wing."""

    S: float
    C_L_alpha: float
    x_ac: float  # over c_bar
    eta: float = 1.0
    downwash_grad: float = 0.0  # d eps / d alpha at the surface (upwash for a canard)
    epsilon_0: float = 0.0
    canard: bool = False

    def effective_slope(self, S_ref):
        """``eta S/S_ref C_L_alpha (1 -+ d eps/d alpha)``."""
        sign = 1 if self.canard else -1
        return self.eta * self.S / S_ref * self.C_L_alpha * (1 + sign * self.downwash_grad)


@dataclass
class LongitudinalStability:
    C_L0: float
    C_L_alpha: float
    C_L_i: float  # per rad of surface incidence
    C_m0: float
    C_m_alpha: float
    C_m_i: float
    x_ac: float  # aircraft neutral point over c_bar
    static_margin: float

    def trim_incidence(self, alpha):
        """Surface incidence for C_m = 0 at ``alpha``."""
        return -(self.C_m0 + self.C_m_alpha * alpha) / self.C_m_i


def longitudinal(S, C_L0_wf, C_L_alpha_wf, x_ac_wf, C_m_ac, x_cg, surface: Surface, ac_shift=0.0):
    """Longitudinal coefficients of wing-body + one horizontal surface (AE 550 final exam questions 12-20).

    The zero-alpha downwash term follows Roskam (tail moment reduced by
    downwash). The part 2 exam script added it with the opposite sign and
    used ``C_l_beta_h`` in place of the tail lift slope in ``C_L0``.
    """
    k = surface.effective_slope(S)
    q = surface.eta * surface.S / S * surface.C_L_alpha
    x_wf = x_ac_wf + ac_shift
    sign = 1 if surface.canard else -1  # downwash lowers a tail's angle of attack, upwash raises a canard's
    C_L0 = C_L0_wf + sign * q * surface.epsilon_0
    C_m0 = C_m_ac + C_L0_wf * (x_cg - x_wf) + sign * q * surface.epsilon_0 * (x_cg - surface.x_ac)
    C_m_alpha = C_L_alpha_wf * (x_cg - x_wf) + k * (x_cg - surface.x_ac)
    x_ac = (C_L_alpha_wf * x_wf + k * surface.x_ac) / (C_L_alpha_wf + k)
    return LongitudinalStability(
        C_L0=C_L0,
        C_L_alpha=C_L_alpha_wf + k,
        C_L_i=q,
        C_m0=C_m0,
        C_m_alpha=C_m_alpha,
        C_m_i=q * (x_cg - surface.x_ac),
        x_ac=x_ac,
        static_margin=x_ac - x_cg,
    )


@dataclass
class LateralStability:
    C_y_beta: float
    C_l_beta: float
    C_n_beta: float
    C_y_dr: float
    C_l_dr: float
    C_n_dr: float

    def rudder_for_sideslip(self, beta):
        """Rudder to zero the yawing moment at sideslip ``beta``."""
        return -self.C_n_beta * beta / self.C_n_dr

    def aileron_for_sideslip(self, beta, C_l_da):
        """Aileron to zero the rolling moment at sideslip ``beta``."""
        return -self.C_l_beta * beta / C_l_da


def lateral(S, b, wf, vertical, horizontal=None):
    """Lateral-directional build-up (AE 550 final exam question 21, part 2 question 8).

    ``wf`` is ``(C_y_beta, C_l_beta, C_n_beta)`` of the wing-body.
    ``vertical`` is a dict with ``S, C_L_alpha, eta, sidewash_grad, x, z, tau_r``
    (``x``, ``z``: fin a.c. aft of and above the c.g., same units as ``b``).
    ``horizontal`` optionally adds ``(S, b, eta, C_y_beta, C_l_beta, C_n_beta)``
    of a tail or canard, scaled by its area and span.
    """
    C_y, C_l, C_n = wf
    if horizontal is not None:
        S_h, b_h, eta_h, cy, cl, cn = horizontal
        C_y += eta_h * S_h / S * cy
        C_l += eta_h * S_h / S * b_h / b * cl
        C_n += eta_h * S_h / S * b_h / b * cn
    v = vertical
    side = v["eta"] * v["S"] / S * v["C_L_alpha"]
    fin = -side * (1 - v["sidewash_grad"])
    rudder = side * v["tau_r"]
    return LateralStability(
        C_y_beta=C_y + fin,
        C_l_beta=C_l + fin * v["z"] / b,
        C_n_beta=C_n - fin * v["x"] / b,
        C_y_dr=rudder,
        C_l_dr=rudder * v["z"] / b,
        C_n_dr=-rudder * v["x"] / b,
    )


def thrust_required(q_bar, S, C_L, C_D, alpha, W, theta):
    """Thrust along body x for steady flight: ``q S (C_D cos a - C_L sin a) + W sin theta``.

    The AE 550 exam scripts used ``+ C_L sin a`` and ``- W sin theta``; with
    lift and drag in wind axes and theta positive nose-up these signs are
    the physical ones.
    """
    return q_bar * S * (C_D * np.cos(alpha) - C_L * np.sin(alpha)) + W * np.sin(theta)
