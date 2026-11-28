"""Orthotropic ply (lamina) used as the building block of a :class:`Laminate`."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

FAILURE_CRITERIA = ("Max_Strain", "Max_Stress")

# Order of the five failure modes in allowables and failure ratios.
_MODE_NAMES = {
    1: "axial tensile",
    2: "transverse tensile",
    3: "axial compressive",
    4: "transverse compressive",
    5: "shear",
}


@dataclass
class PlyFailure:
    """Details of one failure mode reached in one ply."""

    ply_number: int  # 1-based, counted from the top of the laminate
    failure_mode: int  # 1 axial tens., 2 transverse tens., 3 axial comp., 4 transverse comp., 5 shear
    failure_criterion: str  # "Max_Strain" or "Max_Stress"
    failure_load: np.ndarray  # running loads [Nx, Ny, Nxy, Mx, My, Mxy] at which it fails

    def __str__(self) -> str:
        crit = "stress" if self.failure_criterion == "Max_Stress" else "strain"
        loads = ", ".join(f"{v:.3g}" for v in np.ravel(self.failure_load))
        return (
            f"Ply {self.ply_number} has failed in max {_MODE_NAMES[self.failure_mode]} "
            f"{crit} at load [{loads}]"
        )


class Ply:
    """Orthotropic ply for classical lamination theory.

    Parameters
    ----------
    E_1, E_2 : Young's moduli along and across the fibres (Pa).
    G_12 : in-plane shear modulus (Pa).
    nu_12 : major Poisson's ratio.
    stress_allowables : ``[F_1t, F_2t, F_1c, F_2c, tau_12]`` (Pa), optional.
    strain_allowables : ``[eps_1t, eps_2t, eps_1c, eps_2c, gamma_12]``, optional.

    Compressive allowables should be given as negative numbers: a mode
    fails when ``allowable / actual`` lies between 0 and 1.
    """

    def __init__(self, E_1, E_2, G_12, nu_12, stress_allowables=None, strain_allowables=None):
        self.E_1, self.E_2, self.G_12, self.nu_12 = E_1, E_2, G_12, nu_12

        self.failure_stresses_defined = stress_allowables is not None
        self.failure_strains_defined = strain_allowables is not None
        self.stress_allowables = np.asarray(
            stress_allowables if stress_allowables is not None else np.zeros(5), dtype=float
        )
        self.strain_allowables = np.asarray(
            strain_allowables if strain_allowables is not None else np.zeros(5), dtype=float
        )

        # Reduced stiffnesses
        nu_21 = nu_12 * E_2 / E_1
        Q_11 = E_1 / (1 - nu_12 * nu_21)
        Q_22 = E_2 / (1 - nu_12 * nu_21)
        Q_12 = nu_12 * E_2 / (1 - nu_12 * nu_21)
        Q_66 = G_12
        self.Q = np.array([[Q_11, Q_12, 0.0], [Q_12, Q_22, 0.0], [0.0, 0.0, Q_66]])

        # Invariants, used to rotate Q into laminate axes
        self.U_1 = (3 * Q_11 + 3 * Q_22 + 2 * Q_12 + 4 * Q_66) / 8
        self.U_2 = (Q_11 - Q_22) / 2
        self.U_3 = (Q_11 + Q_22 - 2 * Q_12 - 4 * Q_66) / 8
        self.U_4 = (Q_11 + Q_22 + 6 * Q_12 - 4 * Q_66) / 8
        self.U_5 = (Q_11 + Q_22 - 2 * Q_12 + 4 * Q_66) / 8

    # Named accessors for the allowables
    F_1_t = property(lambda self: self.stress_allowables[0])
    F_2_t = property(lambda self: self.stress_allowables[1])
    F_1_c = property(lambda self: self.stress_allowables[2])
    F_2_c = property(lambda self: self.stress_allowables[3])
    epsilon_1_t = property(lambda self: self.strain_allowables[0])
    epsilon_2_t = property(lambda self: self.strain_allowables[1])
    epsilon_1_c = property(lambda self: self.strain_allowables[2])
    epsilon_2_c = property(lambda self: self.strain_allowables[3])

    def Qbar(self, theta: float) -> np.ndarray:
        """Transformed reduced stiffness matrix for a ply at angle ``theta`` (rad)."""
        c2, c4 = np.cos(2 * theta), np.cos(4 * theta)
        s2, s4 = np.sin(2 * theta), np.sin(4 * theta)
        Q11 = self.U_1 + self.U_2 * c2 + self.U_3 * c4
        Q12 = self.U_4 - self.U_3 * c4
        Q22 = self.U_1 - self.U_2 * c2 + self.U_3 * c4
        Q16 = 0.5 * self.U_2 * s2 + self.U_3 * s4
        Q26 = 0.5 * self.U_2 * s2 - self.U_3 * s4
        Q66 = self.U_5 - self.U_3 * c4
        return np.array([[Q11, Q12, Q16], [Q12, Q22, Q26], [Q16, Q26, Q66]])

    def test_failure(self, theta, epsilon, running_loads, failure_criterion, ply_number):
        """Check this ply for failure under laminate-axis strains ``epsilon``.

        Returns ``(has_failed, failures, failure_ratios)`` where
        ``failure_ratios`` is ``allowable / actual`` for the five modes.
        """
        if failure_criterion not in FAILURE_CRITERIA:
            raise NotImplementedError(f"{failure_criterion} failure theory not implemented")
        if (failure_criterion == "Max_Strain" and not self.failure_strains_defined) or (
            failure_criterion == "Max_Stress" and not self.failure_stresses_defined
        ):
            raise ValueError(f"Allowables for {failure_criterion} failure theory not defined")

        c, s = np.cos(theta), np.sin(theta)
        T_2 = np.array(
            [
                [c**2, s**2, -s * c],
                [s**2, c**2, s * c],
                [2 * s * c, -2 * s * c, c**2 - s**2],
            ]
        )
        eps_principal = np.linalg.solve(T_2, np.asarray(epsilon, dtype=float))

        if failure_criterion == "Max_Strain":
            actual, allow = eps_principal, self.strain_allowables
        else:
            actual, allow = self.Q @ eps_principal, self.stress_allowables

        with np.errstate(divide="ignore"):
            ratios = np.array(
                [
                    allow[0] / actual[0],
                    allow[1] / actual[1],
                    allow[2] / actual[0],
                    allow[3] / actual[1],
                    abs(allow[4] / actual[2]),
                ]
            )

        failed_modes = (ratios < 1) & (ratios > 0)
        failures = [
            PlyFailure(ply_number, int(mode) + 1, failure_criterion, ratios[mode] * np.asarray(running_loads))
            for mode in np.flatnonzero(failed_modes)
        ]
        return bool(failed_modes.any()), failures, ratios


def _halpin_tsai(modulus_f, modulus_m, V_f, xi):
    """Halpin-Tsai estimate of a ply's transverse or shear modulus."""
    eta = (modulus_f / modulus_m - 1) / (modulus_f / modulus_m + xi)
    return modulus_m * (1 + xi * eta * V_f) / (1 - eta * V_f)


def mechanical_properties(E_f, G_f, nu_f, E_m, G_m, nu_m, V_f, xi_transverse=None, xi_shear=None):
    """Ply properties from fibre and matrix properties (micromechanics).

    Uses the rule of mixtures for ``E_1`` and ``nu_12`` and Halpin-Tsai
    for ``E_2`` and ``G_12``. ``xi_transverse`` defaults to 2 (circular or
    square fibres; use ``2*a/b`` for rectangular), ``xi_shear`` to 1.

    Returns ``(E_1, E_2, nu_12, G_12)``.
    """
    if xi_transverse is None:
        warnings.warn("xi values not provided: using defaults", stacklevel=2)
        xi_transverse = 2.0
    if xi_shear is None:
        xi_shear = 1.0

    E_1 = E_f * V_f + E_m * (1 - V_f)
    nu_12 = nu_f * V_f + nu_m * (1 - V_f)
    E_2 = _halpin_tsai(E_f, E_m, V_f, xi_transverse)
    G_12 = _halpin_tsai(G_f, G_m, V_f, xi_shear)
    return E_1, E_2, nu_12, G_12


# Matches the MATLAB static method Ply.mechanicalProperties
Ply.mechanical_properties = staticmethod(mechanical_properties)
