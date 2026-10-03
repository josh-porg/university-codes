"""Read and write Tecplot ASCII files in BLOCK data packing.

Ported from the lab's ``importGridFile`` / ``importTecASCIIdata`` /
``OutputTecASCIIdata`` MATLAB tools used for the combustor DMD work.

Connectivity arrays are 0-based in Python (Tecplot files are 1-based;
conversion happens on read and write).
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np


def _read_numbers(path: str | Path, header_lines: int) -> np.ndarray:
    with open(path) as f:
        for _ in range(header_lines):
            f.readline()
        return np.array(f.read().split(), dtype=float)


def read_grid(path: str | Path, header_lines: int, ndim: int, n_nodes: int, n_elements: int):
    """Read a Tecplot FE grid file.

    Returns ``(xyz, connectivity)``: node coordinates ``(n_nodes, ndim)``
    and 0-based element-to-node connectivity ``(n_elements, 2**ndim)``.
    """
    values = _read_numbers(path, header_lines)
    n_xyz = ndim * n_nodes
    xyz = values[:n_xyz].reshape(ndim, n_nodes).T
    nodes_per_element = 2**ndim
    conn = values[n_xyz : n_xyz + n_elements * nodes_per_element]
    connectivity = conn.reshape(n_elements, nodes_per_element).astype(int) - 1
    return xyz, connectivity


def node_to_cell(node_values: np.ndarray, connectivity: np.ndarray) -> np.ndarray:
    """Average node values onto cells (mean of each element's nodes)."""
    return np.asarray(node_values)[connectivity].mean(axis=1)


def read_block_data(
    path: str | Path,
    header_lines: int,
    variables: Sequence[int],
    n_nodes: int,
    n_elements: int,
    cell_centered: Sequence[int] = (),
    connectivity: np.ndarray | None = None,
) -> np.ndarray:
    """Read selected variables from a Tecplot BLOCK-packed data file.

    Parameters
    ----------
    variables : 0-based positions of the variables to extract (counting
        only the variables stored in this file, in file order).
    cell_centered : positions of variables stored at cell centres
        (``n_elements`` values); all others are nodal (``n_nodes`` values).
    connectivity : if given, nodal variables are averaged onto cells so
        every returned column has ``n_elements`` rows.

    Returns an array with one column per requested variable.
    """
    values = _read_numbers(path, header_lines)
    cell_centered = set(cell_centered)

    # Start offset of every stored variable's block
    n_vars = max(max(variables), max(cell_centered, default=-1)) + 1
    sizes = [n_elements if v in cell_centered else n_nodes for v in range(n_vars)]
    starts = np.concatenate([[0], np.cumsum(sizes)])

    columns = []
    for v in variables:
        block = values[starts[v] : starts[v] + sizes[v]]
        if v not in cell_centered and connectivity is not None:
            block = node_to_cell(block, connectivity)
        columns.append(block)
    return np.column_stack(columns)


def write_block_data(
    path: str | Path,
    xyz: np.ndarray,
    values: np.ndarray,
    variable_names: Sequence[str],
    connectivity: np.ndarray,
    title: str = "flowfield.dat",
    zone_type: str | None = None,
    values_per_line: int = 4,
    solution_time: float = 0.0,
) -> None:
    """Write a Tecplot FE zone in BLOCK packing.

    ``values`` has one column per variable. Columns with as many rows as
    elements are written cell-centred; complex values are written as
    ``real_<name>`` and ``imag_<name>`` pairs (e.g. for DMD modes).
    """
    xyz = np.asarray(xyz, dtype=float)
    values = np.asarray(values)
    if values.ndim == 1:
        values = values[:, None]
    names = list(variable_names)
    if np.iscomplexobj(values):
        names = [f"{p}_{n}" for n in names for p in ("real", "imag")]
        values = np.column_stack([part for col in values.T for part in (col.real, col.imag)])

    n_nodes, ndim = xyz.shape
    n_elements = connectivity.shape[0]
    zone_type = zone_type or {2: "FEQuadrilateral", 3: "FEBrick"}[ndim]
    axis_names = ["x", "y", "z"][:ndim]

    def block(arr):
        rows = [
            " ".join(f"{v:20.15E}" for v in arr[i : i + values_per_line])
            for i in range(0, len(arr), values_per_line)
        ]
        return "\n".join(rows) + "\n"

    with open(path, "w") as f:
        f.write(f'TITLE = "{title}"\n')
        f.write("VARIABLES = " + "\n".join(f'"{n}"' for n in axis_names + names) + "\n")
        f.write('ZONE T="zone 1"\n')
        f.write(f"STRANDID=0, SOLUTIONTIME={solution_time}\n")
        f.write(f"Nodes={n_nodes}, Elements={n_elements}, ZONETYPE={zone_type}\n")
        f.write("DATAPACKING=BLOCK\n")
        if values.shape[0] == n_elements and n_elements != n_nodes:
            first, last = ndim + 1, ndim + values.shape[1]
            span = f"{first}-{last}" if last > first else f"{first}"
            f.write(f"VARLOCATION=([{span}]=CELLCENTERED)\n")
        f.write("DT=(" + "SINGLE " * (ndim + values.shape[1]) + ")\n")
        for j in range(ndim):
            f.write(block(xyz[:, j]))
        for j in range(values.shape[1]):
            f.write(block(values[:, j].real))
        for row in np.asarray(connectivity) + 1:
            f.write(" ".join(f"{n:6d}" for n in row) + "\n")
