"""Idealised single-cell wing box with ten booms (AE 725 final project).

Ports ``Stresses``, ``BuckStresses``, ``string_Stresses``,
``wingbox_objective_function``, ``stringy_boi`` and the three constraint
functions (``wingbox_constraints``, ``bucklebox_constraints``,
``stringbox_constraints``).

Boom layout (``y`` up, ``z`` aft, origin at boom 10):

* 1 forward upper spar cap, 2-5 top stringers, 6 aft upper spar cap,
* 7 aft lower spar cap, 8-9 bottom stringers, 10 forward lower spar cap.

Design vector ``x = [A1, A2, A3, A4, t1, t2, t3, t4, (S1..S6)]``: spar-cap
areas (fwd upper, aft upper, aft lower, fwd lower), thicknesses of the top
skin, aft web, bottom skin and forward web, and optionally the six stringer
areas (otherwise ``stringer_area`` each). Units in, lbf, psi.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class WingBoxStresses:
    shear: np.ndarray  # max shear stress in top skin, aft web, bottom skin, fwd web
    cap_axial: np.ndarray  # spar caps 1, 6, 7, 10
    stringer_axial: np.ndarray  # booms 2, 3, 4, 5, 8, 9
    shear_buckling: np.ndarray  # critical shear stress of each panel
    compression_buckling: np.ndarray
    bending_buckling: np.ndarray
    panel_axial: np.ndarray  # mean axial stress of the most loaded top skin panel, aft web, bottom skin, fwd web
    shear_flow: np.ndarray  # q after each boom


@dataclass
class WingBox:
    h_aft: float = 3.07  # L2
    width: float = 17.41  # L3
    h_fwd: float = 5.52  # L4
    V_y: float = -1880.0
    V_z: float = 0.0
    M_y: float = 0.0
    M_z: float = 125853.0
    shear_arm: float = 1.93  # moment arm of V_y about boom 10 in the q10 equation
    stringer_area: float = 0.08
    E: float = 10.1e6
    poisson: float = 0.33
    K_shear: float = 5.7
    K_compression: float = 4.0
    K_bending: float = 36.0
    sigma_allow_t: float = 28000.0
    sigma_allow_c: float = -20000.0
    tau_allow: float = 18000.0

    @property
    def top_length(self):
        return np.hypot(self.h_fwd - self.h_aft, self.width)  # L1

    def lengths(self):
        """``[L1, L2, L3, L4]`` multiplying ``t1..t4`` in the area."""
        return np.array([self.top_length, self.h_aft, self.width, self.h_fwd])

    def area(self, x):
        """Total cross-sectional area and its gradient (``wingbox_objective_function``/``stringy_boi``).

        Like the MATLAB, fixed stringers are not counted.
        """
        x = np.asarray(x, float)
        grad = np.concatenate([np.ones(4), self.lengths(), np.ones(x.size - 8)])
        return float(grad @ x), grad

    def boom_coordinates(self):
        L2, L3, L4 = self.h_aft, self.width, self.h_fwd
        y = np.r_[L4 - np.arange(6) * (L4 - L2) / 5, 0, 0, 0, 0]
        z = np.r_[np.arange(6) * L3 / 5, L3, 2 * L3 / 3, L3 / 3, 0]
        return y, z

    def boom_areas(self, x):
        """Lumped boom areas including the effective skin and web areas.

        The MATLAB ``Stresses`` used a fixed 0.05 in for the top skin in boom
        1; ``BuckStresses`` corrected it to ``t1``, which is used here.
        """
        x = np.asarray(x, float)
        A1, A2, A3, A4, t1, t2, t3, t4 = x[:8]
        S = x[8:14] if x.size >= 14 else np.full(6, self.stringer_area)
        L1, L2, L3, L4 = self.lengths()
        top, bot = L1 / 5 * t1, L3 / 3 * t3
        return np.array([
            L4 * t4 / 6 + top / 2 + A1,
            top + S[0], top + S[1], top + S[2], top + S[3],
            top / 2 + L2 * t2 / 6 + A2,
            L2 * t2 / 6 + bot / 2 + A3,
            bot + S[4], bot + S[5],
            L4 * t4 / 6 + bot / 2 + A4,
        ])  # fmt: skip

    def stresses(self, x) -> WingBoxStresses:
        x = np.asarray(x, float)
        t = x[4:8]
        L1, L2, L3, L4 = self.lengths()
        B = self.boom_areas(x)
        y, z = self.boom_coordinates()
        yc = y - B @ y / B.sum()
        zc = z - B @ z / B.sum()
        Iz, Iy, Iyz = B @ yc**2, B @ zc**2, B @ (yc * zc)
        det = Iy * Iz - Iyz**2

        P = B / det * ((Iy * self.V_y - Iyz * self.V_z) * yc + (Iz * self.V_z - Iyz * self.V_y) * zc)
        cs = np.cumsum(P)
        q10 = (-L2 * L3 * cs[5] + L3 / 3 * L4 * (-cs[6] - cs[7] - cs[8]) + self.V_y * self.shear_arm) / (
            L2 * L3 + L3 / 3 * L4 * 3
        )
        q = q10 + cs

        i1 = int(np.argmax(np.abs(q[:5])))
        i3 = 6 + int(np.argmax(np.abs(q[6:9])))
        shear = np.array([q[i1] / t[0], q[5] / t[1], q[i3] / t[2], q[9] / t[3]])

        sigma = (-(self.M_z * Iy + self.M_y * Iyz) * yc + (self.M_y * Iz + self.M_z * Iyz) * zc) / det
        caps = sigma[[0, 5, 6, 9]]
        stringers = sigma[[1, 2, 3, 4, 7, 8]]

        b = 0.6 * np.array([L1 / 5, L2, L3 / 3, L4])
        plate = np.pi**2 * self.E / (12 * (1 - self.poisson**2)) * (t / b) ** 2
        panels = 0.5 * (sigma + np.roll(sigma, -1))  # panel between boom i and i+1 (10 -> 1)
        j1 = int(np.argmax(np.abs(panels[:5])))
        j3 = 6 + int(np.argmax(np.abs(panels[6:9])))
        # The MATLAB indexed the bottom-skin panel without the +6 offset; it was unused.
        panel_axial = np.array([panels[j1], panels[5], panels[j3], panels[9]])
        return WingBoxStresses(shear, caps, stringers, self.K_shear * plate, self.K_compression * plate,
                               self.K_bending * plate, panel_axial, q)

    def _axial(self, s):
        return np.where(s > 0, s / self.sigma_allow_t, s / self.sigma_allow_c) - 1

    def constraints(self, x, buckling=False):
        """Normalised constraints ``g <= 0``: shear, spar-cap axial, (stringer axial when the
        stringers are design variables), (buckling), and ``-x <= 0``.

        ``wingbox_constraints`` wrote the stress limits as ``1 - allow/stress``;
        the later versions' ``stress/allow - 1`` is used throughout (same sign).
        """
        x = np.asarray(x, float)
        s = self.stresses(x)
        g = [np.abs(s.shear) / self.tau_allow - 1, self._axial(s.cap_axial)]
        if x.size >= 14:
            g.append(self._axial(s.stringer_axial))
        if buckling:
            tau, sig = s.shear, s.panel_axial
            g.append(np.array([
                abs(tau[0]) / s.shear_buckling[0] - 1,
                abs(tau[2]) / s.shear_buckling[2] - 1,
                abs(sig[0]) / s.compression_buckling[0] - 1,
                np.hypot(tau[1] / s.shear_buckling[1], sig[1] / s.bending_buckling[1]) - 1,
                np.hypot(tau[3] / s.shear_buckling[3], sig[3] / s.bending_buckling[3]) - 1,
            ]))  # fmt: skip
        g.append(-x)
        return np.concatenate(g)
