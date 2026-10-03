"""International Standard Atmosphere (AE 445, AE 521, AE 722).

``isa`` replaces MATLAB's ``atmosisa`` (used all over the coursework) and the
hand-written ``getStandardAtmosphericValues`` family. SI units, geopotential
altitude in metres.

The AE 445 originals mixed feet-based lapse constants (``1 - 6.875e-6 h``,
``36089 ft``) with altitudes in metres; the layer formulas below are the
proper SI ones, so values above sea level differ from those scripts.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

G0 = 9.80665  # m/s^2
R_AIR = 287.05287  # J/(kg K)
GAMMA_AIR = 1.4
T0 = 288.15  # K
P0 = 101325.0  # Pa
RHO0 = 1.225  # kg/m^3

# (base altitude m, lapse rate K/m) for the 1976 standard atmosphere layers.
_LAYERS = [
    (0.0, -0.0065),
    (11000.0, 0.0),
    (20000.0, 0.001),
    (32000.0, 0.0028),
    (47000.0, 0.0),
    (51000.0, -0.0028),
    (71000.0, -0.002),
    (84852.0, 0.0),
]


def _layer_bases():
    T, p = [T0], [P0]
    for (h0, L), (h1, _) in zip(_LAYERS[:-1], _LAYERS[1:]):
        T1 = T[-1] + L * (h1 - h0)
        if L == 0:
            p1 = p[-1] * np.exp(-G0 * (h1 - h0) / (R_AIR * T[-1]))
        else:
            p1 = p[-1] * (T1 / T[-1]) ** (-G0 / (L * R_AIR))
        T.append(T1)
        p.append(p1)
    return np.array(T), np.array(p)


_T_BASE, _P_BASE = _layer_bases()


def sutherland_viscosity(T):
    """Dynamic viscosity of air (Pa s) from Sutherland's law (``getMuAir``)."""
    T = np.asarray(T, dtype=float)
    return 1.458e-6 * T**1.5 / (T + 110.4)


def speed_of_sound(T, gamma=GAMMA_AIR, R=R_AIR):
    """Speed of sound (m/s) of a perfect gas at temperature ``T`` (K) (``getSpdSoundAir``)."""
    return np.sqrt(gamma * R * np.asarray(T, dtype=float))


@dataclass
class AtmosphereState:
    T: np.ndarray  # K
    a: np.ndarray  # speed of sound, m/s
    p: np.ndarray  # Pa
    rho: np.ndarray  # kg/m^3
    mu: np.ndarray  # Pa s

    @property
    def theta(self):
        """Temperature ratio T/T0."""
        return self.T / T0

    @property
    def delta(self):
        """Pressure ratio p/p0."""
        return self.p / P0

    @property
    def sigma(self):
        """Density ratio rho/rho0."""
        return self.rho / RHO0

    @property
    def nu(self):
        """Kinematic viscosity (m^2/s)."""
        return self.mu / self.rho


def isa(h, delta_T=0.0) -> AtmosphereState:
    """Standard atmosphere at geopotential altitude ``h`` (m), 0 to 84852 m.

    ``delta_T`` (K) gives an off-standard day at the same pressure altitude,
    like the AE 722 ``Nonstandardrho`` / ``DensityRatio`` helpers: pressure is
    unchanged, temperature is offset, density follows the gas law.
    Altitudes outside the model are clipped to its ends.
    """
    h = np.clip(np.asarray(h, dtype=float), 0.0, _LAYERS[-1][0])
    bases = np.array([b for b, _ in _LAYERS])
    i = np.clip(np.searchsorted(bases, h, side="right") - 1, 0, len(_LAYERS) - 2)
    h0 = bases[i]
    L = np.array([lr for _, lr in _LAYERS])[i]
    Tb, pb = _T_BASE[i], _P_BASE[i]
    T = Tb + L * (h - h0)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(
            L == 0,
            pb * np.exp(-G0 * (h - h0) / (R_AIR * Tb)),
            pb * (T / Tb) ** (-G0 / (np.where(L == 0, 1, L) * R_AIR)),
        )
    T = T + delta_T
    rho = p / (R_AIR * T)
    return AtmosphereState(T=T, a=speed_of_sound(T), p=p, rho=rho, mu=sutherland_viscosity(T))


def density_ratio(h, delta_T=0.0):
    """rho/rho_SL at altitude ``h`` (m) on a day ``delta_T`` (K) off standard."""
    return isa(h, delta_T).rho / RHO0


def true_airspeed(V_indicated, rho, rho_sl=RHO0):
    """Convert equivalent (indicated) to true airspeed (``ias2tas``)."""
    return V_indicated * np.sqrt(rho_sl / np.asarray(rho, dtype=float))


def standard_atmosphere_table(h):
    """Rows of ``(h, T, theta, p, delta, rho, sigma, mu, a)`` for each altitude in ``h``.

    The Python form of ``standardAtmosphericTableGenerator``; feed it to
    :func:`unicodes.io.latex.to_latex_table` for the LaTeX version.
    """
    s = isa(h)
    return np.column_stack(
        [np.asarray(h, dtype=float), s.T, s.theta, s.p, s.delta, s.rho, s.sigma, s.mu, s.a]
    )
