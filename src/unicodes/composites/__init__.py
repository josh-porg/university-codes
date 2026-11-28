"""Classical lamination theory (CLT) for fibre-reinforced composites.

Ported from AE 709 Composites (Autumn 2024)::

    from unicodes.composites import Ply, Laminate

    ply = Ply(E_1=181e9, E_2=10.3e9, G_12=7.17e9, nu_12=0.28)
    lam = Laminate([ply], None, t_ply=0.125e-3,
                   ply_angles=np.deg2rad([0, 90]), symmetry=2)
    lam.A, lam.B, lam.D, lam.E_x
"""

from .laminate import Laminate
from .ply import Ply, PlyFailure, mechanical_properties

__all__ = ["Laminate", "Ply", "PlyFailure", "mechanical_properties"]
