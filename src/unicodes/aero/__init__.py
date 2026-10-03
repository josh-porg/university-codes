"""Aerodynamics and aircraft performance.

    from unicodes.aero import Airfoil, PropellerBlade
    from unicodes.aero import soaring
"""

from . import sizing, soaring
from .airfoil import Airfoil
from .propeller import PropellerBlade, PropellerPerformance

__all__ = ["Airfoil", "PropellerBlade", "PropellerPerformance", "sizing", "soaring"]
