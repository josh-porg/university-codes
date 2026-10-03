"""Modal decompositions of snapshot data: DMD, POD and helpers.

Snapshot convention: one column per time step, i.e. shape
``(n_space, n_time)``. Gridded data with time as the *last* axis, e.g.
``(nx, ny, nt)``, can be flattened with :func:`data_to_snapshots`.

    from unicodes.decomposition import DMD, POD

    dmd = DMD(data, rank=10, dt=1e-4)
    dmd.modes, dmd.omega, dmd.reconstruction
"""

from .dmd import DMD
from .linalg import optimal_svht_coef, rsvd, svht_rank
from .pod import POD
from .snapshots import data_to_snapshots, remove_mean, snapshots_to_data

__all__ = [
    "DMD",
    "POD",
    "data_to_snapshots",
    "optimal_svht_coef",
    "remove_mean",
    "rsvd",
    "snapshots_to_data",
    "svht_rank",
]
