"""Load MATLAB ``.mat`` files into plain Python objects."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def load_mat(path: str | Path) -> dict[str, Any]:
    """Load a ``.mat`` file into a dict of NumPy arrays / nested dicts.

    Equivalent to MATLAB's ``S = load(path)``. Handles both classic
    (v5/v7) files via SciPy and v7.3 (HDF5) files via h5py, which is
    installed with ``pip install "unicodes[mat73]"``.

    MATLAB structs become dicts, so ``S.data.alt`` in MATLAB is
    ``S["data"]["alt"]`` here.
    """
    from scipy.io import loadmat

    path = Path(path)
    try:
        raw = loadmat(path, squeeze_me=True, struct_as_record=False)
    except NotImplementedError:
        return _load_mat73(path)
    return {k: _to_python(v) for k, v in raw.items() if not k.startswith("__")}


def _to_python(value: Any) -> Any:
    from scipy.io.matlab import mat_struct

    if isinstance(value, mat_struct):
        return {name: _to_python(getattr(value, name)) for name in value._fieldnames}
    if isinstance(value, np.ndarray) and value.dtype == object:
        return [_to_python(v) for v in value.ravel()]
    return value


def _load_mat73(path: Path) -> dict[str, Any]:
    try:
        import h5py
    except ImportError as exc:
        raise ImportError(
            f"{path} is a MATLAB v7.3 file; install h5py with "
            '`pip install "unicodes[mat73]"` to read it.'
        ) from exc

    def walk(node: Any) -> Any:
        if isinstance(node, h5py.Group):
            return {k: walk(node[k]) for k in node.keys() if not k.startswith("#")}
        # MATLAB stores arrays column-major, so transpose back.
        return np.array(node).T

    with h5py.File(path, "r") as f:
        return {k: walk(f[k]) for k in f.keys() if not k.startswith("#")}
