"""Modal decompositions of snapshot data: DMD, POD and helpers.

Snapshot convention: one column per time step, i.e. shape
``(n_space, n_time)``. Gridded data with time as the *last* axis, e.g.
``(nx, ny, nt)``, can be flattened with :func:`data_to_snapshots`.

    from unicodes.decomposition import DMD, POD

    dmd = DMD(data, rank=10, dt=1e-4)
    dmd.modes, dmd.omega, dmd.reconstruction
"""

from .dmd import DMD, exact_dmd
from .dmdc import DMDcResult, dmdc
from .linalg import cosamp, fast_ica, optimal_svht_coef, rsvd, svht_rank
from .pod import POD
from .sindy import SINDy, library, library_names, lorenz, stls, stridge, van_der_pol
from .snapshots import data_to_snapshots, remove_mean, snapshots_to_data
from .spectral import dmd_spectrum, fft_amplitude, hankel, mode_power
from .variants import MrDMD, dmd_eigenvalues, era, kernel_dmd, mrdmd, optimal_amplitudes, time_delay_stack

__all__ = [
    "DMD",
    "MrDMD",
    "cosamp",
    "dmd_eigenvalues",
    "era",
    "fast_ica",
    "kernel_dmd",
    "mrdmd",
    "optimal_amplitudes",
    "time_delay_stack",
    "DMDcResult",
    "POD",
    "SINDy",
    "dmd_spectrum",
    "dmdc",
    "exact_dmd",
    "fft_amplitude",
    "hankel",
    "library",
    "library_names",
    "lorenz",
    "mode_power",
    "stls",
    "stridge",
    "van_der_pol",
    "data_to_snapshots",
    "optimal_svht_coef",
    "remove_mean",
    "rsvd",
    "snapshots_to_data",
    "svht_rank",
]
