"""Ride-discomfort index for a rigid wing vs a passive-aeroelastic (PAH) wing (``Ride_quality``)."""

import numpy as np

from unicodes.aero.soaring import ride_quality_index
from unicodes.atmosphere import isa
from unicodes.units import FT

rho = float(isa(3000 * FT, 10).rho)
C_L_alpha = np.array([2 * np.pi, 2 * np.pi / 4])  # rigid flat plate, PAH with 4x relief
print("C_ride [rigid, PAH] =", ride_quality_index(rho, 28.2, 500.0, C_L_alpha, -0.248))
