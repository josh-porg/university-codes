"""Thin-airfoil theory for an arbitrary camber line (AE 546 HW 4).

The MATLAB homework integrated the Fourier coefficients symbolically for a
piecewise camber line; here they are integrated numerically, so any
``dz/dx`` function works.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad


@dataclass
class ThinAirfoil:
    camber_slope: callable  # dz/dx as a function of x/c in [0, 1]
    breakpoints: tuple = ()  # x/c locations where the slope is discontinuous

    def _theta_points(self):
        return [np.arccos(1 - 2 * x) for x in self.breakpoints]

    def _integrate(self, f):
        g = lambda th: self.camber_slope((1 - np.cos(th)) / 2) * f(th)  # noqa: E731
        return quad(g, 0, np.pi, points=self._theta_points() or None, limit=200)[0]

    def fourier_coefficient(self, n: int, alpha: float = 0.0) -> float:
        """A_n of the vortex-sheet series (A_0 depends on ``alpha``)."""
        if n == 0:
            return alpha - self._integrate(lambda th: 1.0) / np.pi
        return 2 / np.pi * self._integrate(lambda th: np.cos(n * th))

    @property
    def zero_lift_angle(self) -> float:
        """alpha_L=0 = -(1/pi) int dz/dx (cos(theta) - 1) dtheta."""
        return -self._integrate(lambda th: np.cos(th) - 1) / np.pi

    def c_l(self, alpha):
        return 2 * np.pi * (np.asarray(alpha) - self.zero_lift_angle)

    @property
    def c_m_quarter_chord(self) -> float:
        return np.pi / 4 * (self.fourier_coefficient(2) - self.fourier_coefficient(1))

    def center_of_pressure(self, alpha) -> float:
        """x_cp/c = 1/4 (1 + pi/c_l (A1 - A2))."""
        A1, A2 = self.fourier_coefficient(1), self.fourier_coefficient(2)
        return 0.25 * (1 + np.pi / self.c_l(alpha) * (A1 - A2))


def two_piece_parabolic_camber(x_split, a_fwd, b_fwd, c_fwd, a_aft, b_aft, c_aft):
    """Thin airfoil for ``z/c = a (b + c x - x^2)``-type pieces either side of ``x_split``.

    AE 546 HW 4.6 uses ``0.25 (0.8 x - x^2)`` for x < 0.4 and
    ``0.111 (0.2 + 0.8 x - x^2)`` aft: ``two_piece_parabolic_camber(0.4, 0.25, 0, 0.8, 0.111, 0.2, 0.8)``.
    """

    def slope(x):
        return np.where(x < x_split, a_fwd * (c_fwd - 2 * x), a_aft * (c_aft - 2 * x))

    return ThinAirfoil(slope, (x_split,))
