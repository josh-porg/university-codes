"""Blade-element momentum analysis of the Archytas propeller (``Blade_Element_Momentum_Theory``).

A port of Barrett's spreadsheet; see :class:`unicodes.aero.PropellerBlade`.
"""

import numpy as np

from unicodes.aero import PropellerBlade
from unicodes.units import KT, SLUG_PER_FT3

blade = PropellerBlade(np.deg2rad(49), np.deg2rad(16), np.deg2rad(-40), -0.7, root_chord=0.07, radius=0.9)
perf = blade.analyze(omega=1439 * 2 * np.pi / 60, V=68 * KT, rho=0.002133 * SLUG_PER_FT3)
sigma = blade.n_blades * blade.root_chord * (1 + blade.taper) / 2 / (np.pi * blade.radius)
print(f"C_T = {perf.C_T:.5f}, C_P = {perf.C_P:.5f}, C_T/sigma = {perf.C_T / sigma:.4f}")
print(f"T = {perf.thrust:.1f} N, P_shaft = {perf.shaft_power:.0f} W, eta_p = {perf.efficiency:.3f}")
