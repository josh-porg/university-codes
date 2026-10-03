"""Airfoil lift curve from a few key points (AE 722 sailplane design).

The lift curve is a clamped cubic spline through zero lift, the end of
the linear range, c_l_max and a post-stall point. Extra points along the
linear range keep it straight, and the end slope is chosen so the curve
peaks exactly at ``alpha_c_l_max``. All angles in radians.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass(frozen=True)
class Airfoil:
    name: str
    data_source: str
    c_l_alpha: float  # lift-curve slope (1/rad)
    alpha_0: float  # zero-lift angle (rad)
    c_l_0: float  # lift coefficient at zero angle
    alpha_star: float  # end of the linear range (rad)
    alpha_prime: float  # (rad), stored for the airfoil file
    c_l_star: float  # c_l at alpha_star
    alpha_c_l_max: float  # rad
    c_l_max: float
    alpha_post_stall: float  # rad
    c_l_post_stall: float
    linear_points: int = 2  # extra points forcing linearity (1-5)

    def key_points(self) -> tuple[np.ndarray, np.ndarray]:
        """``(alpha, c_l)`` of the four defining points."""
        alpha = np.array([self.alpha_0, self.alpha_star, self.alpha_c_l_max, self.alpha_post_stall])
        c_l = np.array([0.0, self.c_l_star, self.c_l_max, self.c_l_post_stall])
        return alpha, c_l

    @cached_property
    def spline(self) -> CubicSpline:
        alpha, c_l = self.key_points()
        n = self.linear_points
        frac = np.arange(1, n + 1) / (n + 1)
        x = np.concatenate([[alpha[0]], alpha[0] + frac * (alpha[1] - alpha[0]), alpha[1:]])
        y = np.concatenate([[c_l[0]], c_l[0] + frac * (c_l[1] - c_l[0]), c_l[1:]])
        start_slope = (c_l[1] - c_l[0]) / (alpha[1] - alpha[0])

        def make(end_slope):
            return CubicSpline(x, y, bc_type=((1, start_slope), (1, end_slope)))

        # The slope at alpha_c_l_max is linear in the end slope, so two
        # evaluations give the end slope that puts the peak there exactly
        # (the MATLAB iterated on a scale factor to the same effect).
        d0 = make(0.0)(alpha[2], 1)
        d1 = make(1.0)(alpha[2], 1)
        return make(-d0 / (d1 - d0))

    def c_l(self, alpha):
        """Lift coefficient; zero beyond the post-stall point."""
        alpha = np.asarray(alpha, dtype=float)
        return np.where(alpha <= self.alpha_post_stall, self.spline(alpha), 0.0)

    def c_l_alpha_at(self, alpha):
        """Local lift-curve slope dc_l/dalpha."""
        return self.spline(np.asarray(alpha, dtype=float), 1)

    def write_airfoil_file(self, path: str | Path) -> None:
        """Write the ``[Airfoil]`` text file format used by the design tools."""
        deg = np.rad2deg
        Path(path).write_text(
            "[Airfoil]\n"
            f"Name={self.name}\n"
            f"Data Source={self.data_source}\n"
            f"clalpha [1/rad]={self.c_l_alpha:.14f}\n"
            f"alpha_o [deg]={deg(self.alpha_0):.14f}\n"
            f"c_l_o [-]={self.c_l_0:.14f}\n"
            f'alpha" [deg]={deg(self.alpha_prime):.14f}\n'
            f"alpha* [deg]={deg(self.alpha_star):.14f}\n"
            f"c_l* [-]={self.c_l_star:.14f}\n"
            f"alpha_c_l_max [deg]={deg(self.alpha_c_l_max):.14f}\n"
            f"c_l_max [-]={self.c_l_max:.14f}\n"
            f"alpha_post_S [deg]={deg(self.alpha_post_stall):.14f}\n"
            f"c_l_post_S [-]={self.c_l_post_stall:.14f}\n"
        )
