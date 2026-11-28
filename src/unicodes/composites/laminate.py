"""Laminate stiffness (ABD) and first-ply failure via classical lamination theory.

Assumes thin laminates in plane stress where a line straight and normal
to the mid-plane stays straight, normal and the same length (Kirchhoff).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .ply import Ply, PlyFailure

_SYMMETRIES = {
    0: "asymmetric",
    1: "odd symmetric",
    2: "even symmetric",
    -1: "odd anti-symmetric",
    -2: "even anti-symmetric",
}


class Laminate:
    """A stack of plies analysed with classical lamination theory.

    Parameters
    ----------
    plies : the distinct :class:`Ply` materials used.
    ply_material : for each layer (top to bottom) the index into ``plies``
        of its material. ``None`` means every layer uses ``plies[0]``.
        Unlike MATLAB, indices are 0-based.
    t_ply : ply thickness (m), a scalar or one value per layer.
    ply_angles : angle of each layer (rad), top to bottom, within ±pi/2.
    symmetry : how to mirror the given layers: 0 none, 1 odd symmetric,
        2 even symmetric, -1 odd anti-symmetric, -2 even anti-symmetric.
        Odd symmetry mirrors about the last given ply.
    """

    def __init__(
        self,
        plies: Sequence[Ply],
        ply_material: Sequence[int] | None,
        t_ply,
        ply_angles,
        symmetry: int = 0,
    ):
        ply_angles = np.ravel(np.asarray(ply_angles, dtype=float))
        if np.any(np.abs(ply_angles) > np.pi / 2):
            raise ValueError(
                f"Ply angles {ply_angles[np.abs(ply_angles) > np.pi / 2]} exceed ±pi/2 rad; "
                "were they given in degrees?"
            )
        if symmetry not in _SYMMETRIES:
            raise ValueError(f"symmetry {symmetry} not valid; must be one of {_SYMMETRIES}")
        if abs(symmetry) == 1 and ply_angles[-1] not in (0.0, np.pi / 2):
            if symmetry == -1:
                raise ValueError("Odd anti-symmetry about an angle ply is impossible")
            import warnings

            warnings.warn("Odd symmetry about an angle ply may give misleading failure results", stacklevel=2)

        n = len(ply_angles)
        material = np.zeros(n, dtype=int) if ply_material is None else np.ravel(ply_material).astype(int)
        if len(material) != n:
            raise ValueError(f"ply_material has {len(material)} entries but there are {n} angles")
        t = np.full(n, float(t_ply)) if np.isscalar(t_ply) else np.ravel(np.asarray(t_ply, dtype=float))
        if len(t) != n:
            raise ValueError(f"t_ply has {len(t)} entries but there are {n} angles")

        # Mirror the given half of the laminate
        if symmetry == 0:
            mirror = slice(0, 0)
        elif abs(symmetry) == 1:
            mirror = slice(0, n - 1)
        else:
            mirror = slice(0, n)
        mirrored_angles = ply_angles[mirror][::-1]
        if symmetry < 0:
            mirrored_angles = -mirrored_angles
            mirrored_angles[np.isclose(mirrored_angles, -np.pi / 2)] = np.pi / 2

        self.plies = list(plies)
        self.symmetry = symmetry
        self.ply_angles = np.concatenate([ply_angles, mirrored_angles])
        self.t_ply = np.concatenate([t, t[mirror][::-1]])
        self.ply_material = np.concatenate([material, material[mirror][::-1]])
        self.n_plies = len(self.ply_angles)

        # Ply boundary positions, measured from the mid-plane, top surface most negative
        z = np.concatenate([[0.0], np.cumsum(self.t_ply)])
        self.t_laminate = z[-1]
        self.z = z - self.t_laminate / 2

        A = np.zeros((3, 3))
        B = np.zeros((3, 3))
        D = np.zeros((3, 3))
        for k in range(self.n_plies):
            Qbar = self.plies[self.ply_material[k]].Qbar(self.ply_angles[k])
            z0, z1 = self.z[k], self.z[k + 1]
            A += Qbar * (z1 - z0)
            B += 0.5 * Qbar * (z1**2 - z0**2)
            D += (1 / 3) * Qbar * (z1**3 - z0**3)

        self.laminate_stiffness = np.block([[A, B], [B, D]])
        self.laminate_compliance = np.linalg.inv(self.laminate_stiffness)

        a = self.laminate_compliance
        self.E_x = 1 / (self.t_laminate * a[0, 0])
        self.E_y = 1 / (self.t_laminate * a[1, 1])
        self.G_xy = 1 / (self.t_laminate * a[2, 2])
        self.nu_xy = -a[0, 1] / a[0, 0]

    @property
    def A(self) -> np.ndarray:
        """Extensional stiffness matrix."""
        return self.laminate_stiffness[:3, :3]

    @property
    def B(self) -> np.ndarray:
        """Extension-bending coupling matrix."""
        return self.laminate_stiffness[:3, 3:]

    @property
    def D(self) -> np.ndarray:
        """Bending stiffness matrix."""
        return self.laminate_stiffness[3:, 3:]

    def mid_plane_response(self, running_loads) -> tuple[np.ndarray, np.ndarray]:
        """Mid-plane strains ``epsilon_0`` and curvatures ``kappa`` under
        running loads ``[Nx, Ny, Nxy, Mx, My, Mxy]``."""
        response = self.laminate_compliance @ np.asarray(running_loads, dtype=float)
        return response[:3], response[3:]

    def load(self, running_loads, failure_criterion: str = "Max_Strain"):
        """Apply running loads and check every ply for failure.

        Each ply is checked at its top and bottom surface.

        Returns
        -------
        failed_plies : bool array, one entry per ply.
        failures : list of :class:`PlyFailure`.
        failure_ratios : array (n_plies, 2 surfaces, 5 modes) of allowable/actual.
        """
        running_loads = np.asarray(running_loads, dtype=float)
        epsilon_0, kappa = self.mid_plane_response(running_loads)

        failed_plies = np.zeros(self.n_plies, dtype=bool)
        failures: list[PlyFailure] = []
        ratios = np.zeros((self.n_plies, 2, 5))
        for surface in range(2):  # 0 = top of each ply, 1 = bottom
            for k in range(self.n_plies):
                eps = epsilon_0 + self.z[k + surface] * kappa
                ply = self.plies[self.ply_material[k]]
                has_failed, ply_failures, r = ply.test_failure(
                    self.ply_angles[k], eps, running_loads, failure_criterion, k + 1
                )
                failed_plies[k] |= has_failed
                failures.extend(ply_failures)
                ratios[k, surface] = r
        return failed_plies, failures, ratios
