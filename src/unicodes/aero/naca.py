"""NACA 4- and 5-digit airfoil coordinates (AE 546 ``naca5gen``).

``naca5gen`` was Divahar Jayaraman's generator; this is a port of the same
equations (standard 230-series camber, thickness with closed or open TE).
Coordinates are normalised by chord.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

_THICKNESS = (0.2969, -0.1260, -0.3516, 0.2843)


@dataclass
class AirfoilCoordinates:
    name: str
    x: np.ndarray  # TE -> upper -> LE -> lower -> TE
    z: np.ndarray
    x_camber: np.ndarray
    z_camber: np.ndarray
    le_radius: float

    @property
    def le_index(self) -> int:
        return int(np.argmin(self.x))

    @property
    def upper(self):
        """Upper surface ``(x, z)`` from trailing edge to leading edge."""
        i = self.le_index
        return self.x[: i + 1], self.z[: i + 1]

    @property
    def lower(self):
        """Lower surface ``(x, z)`` from leading edge to trailing edge."""
        i = self.le_index
        return self.x[i:], self.z[i:]

    def write_dat(self, path: str | Path) -> None:
        """Selig-style ``.dat`` file: name line then ``x z`` rows."""
        rows = "\n".join(f"{a:.6f} {b:.6f}" for a, b in zip(self.x, self.z))
        Path(path).write_text(f"{self.name}\n{rows}\n")


def _stations(n: int, cosine: bool) -> np.ndarray:
    if cosine:
        return 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    return np.linspace(0, 1, n + 1)


def _thickness(x, t, finite_te):
    a4 = -0.1015 if finite_te else -0.1036
    a0, a1, a2, a3 = _THICKNESS
    return t / 0.2 * (a0 * np.sqrt(x) + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4)


def _assemble(name, x, yt, zc, slope, t):
    theta = np.arctan(slope)
    xu, zu = x - yt * np.sin(theta), zc + yt * np.cos(theta)
    xl, zl = x + yt * np.sin(theta), zc - yt * np.cos(theta)
    return AirfoilCoordinates(
        name=name,
        x=np.concatenate([xu[::-1], xl[1:]]),
        z=np.concatenate([zu[::-1], zl[1:]]),
        x_camber=x,
        z_camber=zc,
        le_radius=0.5 * (0.2969 * t / 0.2) ** 2,
    )


def naca4(designation: str, n: int = 50, cosine_spacing: bool = True, finite_te: bool = False):
    """NACA 4-digit section with ``n`` panels per surface."""
    m, p, t = int(designation[0]) / 100, int(designation[1]) / 10, int(designation[2:]) / 100
    x = _stations(n, cosine_spacing)
    yt = _thickness(x, t, finite_te)
    if m == 0 or p == 0:
        zc, slope = np.zeros_like(x), np.zeros_like(x)
    else:
        fwd = x <= p
        zc = np.where(fwd, m / p**2 * (2 * p * x - x**2), m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x**2))
        slope = np.where(fwd, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    return _assemble(f"NACA {designation}", x, yt, zc, slope, t)


def naca5(designation: str, n: int = 50, cosine_spacing: bool = True, finite_te: bool = False):
    """NACA 5-digit (standard, non-reflexed) section with ``n`` panels per surface (``naca5gen``)."""
    cl_design = int(designation[0]) * 1.5 / 10
    p = 0.5 * int(designation[1:3]) / 100
    t = int(designation[3:]) / 100
    x = _stations(n, cosine_spacing)
    yt = _thickness(x, t, finite_te)
    if p == 0:
        zc, slope = np.zeros_like(x), np.zeros_like(x)
    else:
        P = [0.05, 0.1, 0.15, 0.2, 0.25]
        M = [0.0580, 0.1260, 0.2025, 0.2900, 0.3910]
        K = [361.4, 51.64, 15.957, 6.643, 3.230]
        m = float(CubicSpline(P, M, bc_type="not-a-knot")(p))
        k1 = float(CubicSpline(M, K, bc_type="not-a-knot")(m))
        fwd = x <= p
        scale = cl_design / 0.3
        zc = scale * np.where(fwd, k1 / 6 * (x**3 - 3 * m * x**2 + m**2 * (3 - m) * x), k1 / 6 * m**3 * (1 - x))
        # The original left the slope unscaled and gave the aft slope a + sign;
        # both only matter away from the 230 series.
        slope = scale * np.where(fwd, k1 / 6 * (3 * x**2 - 6 * m * x + m**2 * (3 - m)), -k1 / 6 * m**3)
    return _assemble(f"NACA {designation}", x, yt, zc, slope, t)


def naca(designation: str, **kwargs):
    """Dispatch on the number of digits."""
    return {4: naca4, 5: naca5}[len(designation)](designation, **kwargs)
