"""Aerodynamics and aircraft performance.

    from unicodes.aero import Airfoil, PropellerBlade
    from unicodes.aero import soaring
"""

from . import soaring
from .airfoil import Airfoil
from .propeller import PropellerBlade, PropellerPerformance

__all__ = ["Airfoil", "PropellerBlade", "PropellerPerformance", "soaring"]
