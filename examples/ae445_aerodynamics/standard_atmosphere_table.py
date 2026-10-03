"""Standard atmosphere tables in SI and British units (``standardAtmosphericTableGenerator``).

Writes ``StandardAtmosphereMTable.tex`` and ``StandardAtmosphereBGTable.tex``.
(The MATLAB BG table had a bug: it overwrote pressure with temperature.)
"""

import numpy as np

from unicodes.atmosphere import standard_atmosphere_table
from unicodes.io import to_latex_table
from unicodes.units import FT, LBF_S_PER_FT2, PSF, RANKINE, SLUG_PER_FT3

HEADERS = ["h", "T", r"$\theta$", "p", r"$\delta$", r"$\rho$", r"$\sigma$", r"$\mu$", "a"]

si = standard_atmosphere_table(np.arange(0, 20001, 500))
print(to_latex_table(si, HEADERS, path="StandardAtmosphereMTable.tex"))

bg = standard_atmosphere_table(np.arange(0, 20000, 1000 * FT))
bg[:, 0] /= FT
bg[:, 1] /= RANKINE
bg[:, 3] /= PSF
bg[:, 5] /= SLUG_PER_FT3
bg[:, 7] /= LBF_S_PER_FT2
bg[:, 8] /= FT
print(to_latex_table(bg, HEADERS, path="StandardAtmosphereBGTable.tex"))
