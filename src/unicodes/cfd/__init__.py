"""Compressible-flow CFD building blocks (AE 746).

    from unicodes.cfd import StructuredMesh2D, solve_steady, euler
"""

from . import advection, euler
from .mesh import StructuredMesh2D
from .navier_stokes import Viscosity
from .solver2d import SteadyResult, residual, solve_steady, solve_unsteady

__all__ = ["StructuredMesh2D", "SteadyResult", "advection", "euler", "Viscosity", "residual", "solve_steady", "solve_unsteady"]
