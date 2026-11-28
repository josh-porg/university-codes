"""Run Mark Drela's XFOIL and read its polar (AE 722 ``xfoil.m``, ``xfoilCl.m``, ``Lift_mesh_Xfoil``).

The MATLAB wrapper (Rafael Oliveira's, from File Exchange) wrote an input
deck, ran ``xfoil.exe`` and parsed the accumulated polar. This does the same
with :mod:`subprocess`; XFOIL must be installed (``xfoil`` on the PATH, or
pass ``executable``). Only the polar is returned.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Polar:
    name: str
    Re: float
    Ncrit: float
    alpha: np.ndarray  # deg
    CL: np.ndarray
    CD: np.ndarray
    CDp: np.ndarray
    Cm: np.ndarray
    top_xtr: np.ndarray
    bot_xtr: np.ndarray


def parse_polar(text: str) -> Polar:
    """Parse an XFOIL ``pacc`` polar file."""
    name = re.search(r"Calculated polar for:\s*(.*)", text)
    re_match = re.search(r"Re\s*=\s*([\d.]+)\s*e\s*([+-]?\d+)", text)
    ncrit = re.search(r"Ncrit\s*=\s*([\d.]+)", text)
    rows = []
    started = False
    for line in text.splitlines():
        if line.strip().startswith("------"):
            started = True
            continue
        if started and line.strip():
            rows.append([float(v) for v in line.split()[:7]])
    data = np.array(rows, dtype=float).reshape(-1, 7)
    return Polar(
        name=name.group(1).strip() if name else "",
        Re=float(re_match.group(1)) * 10 ** int(re_match.group(2)) if re_match else 0.0,
        Ncrit=float(ncrit.group(1)) if ncrit else np.nan,
        alpha=data[:, 0], CL=data[:, 1], CD=data[:, 2], CDp=data[:, 3], Cm=data[:, 4],
        top_xtr=data[:, 5], bot_xtr=data[:, 6],
    )


def run_xfoil(airfoil, alpha, Re=1e6, mach=0.0, commands=(), executable="xfoil", iterations=None, timeout=300):
    """Run XFOIL over the angles ``alpha`` (deg) and return the :class:`Polar`.

    ``airfoil`` is a NACA designation (``"NACA2412"``), a coordinate file
    path, or an ``(n, 2)`` array running TE -> upper -> LE -> lower -> TE.
    ``commands`` are extra XFOIL command strings run before ``oper`` (e.g.
    ``"ppar n 160"``, ``"gdes flap 0.75 0 5 exec"``); spaces or slashes
    separate menu levels like in the MATLAB wrapper.
    """
    exe = shutil.which(executable) or executable
    alpha = np.atleast_1d(alpha)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        lines = []
        if isinstance(airfoil, str) and re.fullmatch(r"(?i)NACA\s*\d{4,5}", airfoil.strip()):
            lines.append(f"naca {airfoil.strip()[4:].strip()}")
        elif isinstance(airfoil, (str, Path)):
            lines.append(f"load {Path(airfoil).resolve()}")
        else:
            coords = tmp / "foil.dat"
            coords.write_text("foil\n" + "\n".join(f"{x:9.5f} {y:9.5f}" for x, y in np.asarray(airfoil)))
            lines.append(f"load {coords}")
        for c in commands:
            lines += re.split(r"[ \\/]+", c) + [""]
        lines += ["", "oper"]
        if iterations:
            lines.append(f"iter {iterations}")
        lines += [f"re {Re:g}", f"mach {mach:g}"]
        if Re > 0:
            lines.append("visc")
        lines += ["pacc", "polar.dat", "", *[f"alfa {a:g}" for a in alpha], "pacc", "", "quit"]
        subprocess.run([exe], input="\n".join(lines) + "\n", text=True, cwd=tmp, capture_output=True, timeout=timeout, check=True)
        polar = parse_polar((tmp / "polar.dat").read_text())
    if len(polar.alpha) != len(alpha):
        import warnings

        warnings.warn("XFOIL did not converge at every angle; try more iterations", stacklevel=2)
    return polar


def lift_mesh(cambers_percent, alpha_max=45.0, n_alpha=191, Re=5e5, run=run_xfoil, **kwargs):
    """C_l(alpha, camber) mesh from NACA x412 sections, mirrored for negative camber (``Lift_mesh_Xfoil``).

    Runs each section from 0 deg up and from 0 deg down (for convergence),
    returns ``(alpha_rad, camber, c_l)`` arrays of shape ``(2 n_alpha - 1, 2 n_camber - 1)``.
    """
    up = np.linspace(0, alpha_max, n_alpha)
    down = np.linspace(0, -alpha_max, n_alpha)
    alpha = np.concatenate([down[:0:-1], up])
    cols = []
    for m in cambers_percent:
        name = f"NACA{int(m)}412"
        p_up, p_dn = run(name, up, Re, **kwargs), run(name, down, Re, **kwargs)
        a = np.concatenate([p_dn.alpha[:0:-1], p_up.alpha])
        cl = np.concatenate([p_dn.CL[:0:-1], p_up.CL])
        cols.append(np.interp(alpha, a, cl, left=np.nan, right=np.nan))
    cl_pos = np.column_stack(cols)
    camber_pos = np.asarray(cambers_percent, dtype=float) / 100
    cl = np.hstack([-cl_pos[::-1, :0:-1], cl_pos])
    camber = np.concatenate([-camber_pos[:0:-1], camber_pos])
    A, C = np.meshgrid(np.deg2rad(alpha), camber, indexing="ij")
    return A, C, cl
