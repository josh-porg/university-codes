"""Lift and pitching-moment surfaces over angle of attack and camber for passive
aero-compliant hinged (PAH) wings (AE 722 ``PAH_Lift_Surface_Gradients``,
``PAH_Moment_Surface_Gradients``).

A PAH wing changes camber with angle of attack. The section characteristics
between a slightly cambered and a highly cambered airfoil are blended
linearly into a surface ``c(alpha, camber)``; a camber schedule
``camber(alpha)`` cuts a curve out of it, which is the PAH wing's effective
lift or moment curve.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import griddata

from .airfoil import Airfoil
from .wing import cm_alpha_tail, cm_alpha_wing_fuselage, cm_tail_incidence


def _resample(airfoil: Airfoil, counts):
    """Points along the lift curve, ``counts`` per region (linear, pre-stall, post-stall)."""
    a, _ = airfoil.key_points()
    alphas = [a[:1]]
    for i, n in enumerate(counts):
        alphas.append(np.linspace(a[i], a[i + 1], n + 1)[1:])
    alpha = np.concatenate(alphas)
    return alpha, airfoil.spline(alpha)


def mesh_gradient(values, alpha, camber):
    """Partial derivatives on the (alpha rows, camber columns) mesh (``MeshGradient``).

    Average of the two one-sided differences inside, one-sided at the edges.
    """
    def along(v, x, axis):
        dv, dx = np.diff(v, axis=axis), np.diff(x, axis=axis)
        s = dv / dx
        out = np.empty_like(v, dtype=float)
        if axis == 0:
            out[1:-1] = (s[1:] + s[:-1]) / 2
            out[0], out[-1] = s[0], s[-1]
        else:
            out[:, 1:-1] = (s[:, 1:] + s[:, :-1]) / 2
            out[:, 0], out[:, -1] = s[:, 0], s[:, -1]
        return out

    return along(values, alpha, 0), along(values, camber, 1)


@dataclass
class CamberSurface:
    alpha: np.ndarray  # (n_alpha, n_camber), rad
    camber: np.ndarray
    values: dict  # name -> (n_alpha, n_camber) array, e.g. {"c_l": ..., "c_m": ...}

    @classmethod
    def from_airfoils(cls, neutral: Airfoil, cambered: Airfoil, camber_neutral, camber_max,
                      counts=(10, 5, 2), n_camber=10, extra=None):
        """Blend two airfoils' lift curves; ``extra`` maps names to ``f(airfoil, alpha)`` evaluated on each."""
        a_n, cl_n = _resample(neutral, counts)
        a_c, cl_c = _resample(cambered, counts)
        w = np.linspace(0, 1, n_camber)
        alpha = a_n[:, None] + (a_c - a_n)[:, None] * w
        camber = np.broadcast_to(camber_neutral + (camber_max - camber_neutral) * w, alpha.shape).copy()
        values = {"c_l": cl_n[:, None] + (cl_c - cl_n)[:, None] * w}
        for name, f in (extra or {}).items():
            v_n, v_c = np.asarray(f(neutral, a_n)), np.asarray(f(cambered, a_c))
            values[name] = v_n[:, None] + (v_c - v_n)[:, None] * w
        return cls(alpha, camber, values)

    def __call__(self, name, alpha, camber):
        """Interpolate ``name`` at scattered ``(alpha, camber)`` points (``griddata``, linear)."""
        pts = np.column_stack([self.alpha.ravel(), self.camber.ravel()])
        return griddata(pts, self.values[name].ravel(), (np.asarray(alpha), np.asarray(camber)), method="linear")

    def slice(self, name, alpha, camber_schedule):
        """Effective curve ``name(alpha, camber(alpha))`` for a camber schedule callable."""
        alpha = np.asarray(alpha, dtype=float)
        return self(name, alpha, camber_schedule(alpha))

    def gradients(self, name):
        """``(d/dalpha, d/dcamber)`` of ``name`` on the mesh."""
        return mesh_gradient(self.values[name], self.alpha, self.camber)

    def slope_along_schedule(self, name, camber_slope):
        """Total derivative ``d name/d alpha = d/dalpha + d/dcamber * dcamber/dalpha`` on the mesh."""
        da, dc = self.gradients(name)
        return da + dc * camber_slope


def linear_camber_schedule(alpha_c, camber_c, slope):
    """``camber(alpha) = camber_c + slope (alpha - alpha_c)``: the planar slices in the MATLAB plots."""
    return lambda alpha: camber_c + slope * (np.asarray(alpha) - alpha_c)


def pitching_moment(alpha, airfoil: Airfoil, i_h, C_m0_wf, C_m0_h, C_L_alpha_h, x_cg, x_ac_wf, x_ac_h,
                    S_h, S, eta_h, downwash_grad):
    """Aircraft C_m by integrating C_m_alpha from 0 to ``alpha`` (``PitchingMomentCoeffcient``).

    ``eta_h`` is a callable tail dynamic-pressure ratio ``eta_h(alpha)``
    (the course fitted it to wind-tunnel blanking data and added 0.2).
    """
    def dcm(a):
        e = eta_h(a)
        return (cm_alpha_wing_fuselage(float(airfoil.c_l_alpha_at(a)), x_cg, x_ac_wf)
                + cm_alpha_tail(S_h, S, e, C_L_alpha_h, downwash_grad, x_cg, x_ac_h)
                + cm_tail_incidence(S_h, S, e, C_L_alpha_h, x_cg, x_ac_h) * i_h)

    return np.array([C_m0_wf + C_m0_h + quad(dcm, 0, a, limit=200)[0] for a in np.atleast_1d(alpha)])
