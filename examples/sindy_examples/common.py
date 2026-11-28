"""Shared pieces of the SINDy example ports (``poolData``, ``poolDataLIST``, ``sparsifyDynamics``,
``sparseGalerkin``, ``color_line3``) on top of :mod:`unicodes.decomposition.sindy`.

Code by S. L. Brunton for Brunton, Proctor & Kutz, "Discovering governing
equations from data by sparse identification of nonlinear dynamical
systems", PNAS 113 (2016).
"""

import numpy as np
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.integrate import solve_ivp

from unicodes.decomposition import library, library_names, stls


def terms(polyorder, usesine=False):
    """Library terms of ``poolData(..., polyorder, usesine)``."""
    t = ["const", "linear", "poly2", "poly3", "poly4", "poly5"][: polyorder + 1]
    return tuple(t + (["sin"] if usesine else []))


def pool_data(X, polyorder, usesine=False):
    return library(np.atleast_2d(X), terms(polyorder, usesine))


def print_model(Xi, names, polyorder, usesine=False):
    """The coefficient table of ``poolDataLIST``."""
    labels = library_names(names, terms(polyorder, usesine))
    print(f"{'':10s}" + "".join(f"{'d' + n + '/dt':>14s}" for n in names))
    for lab, row in zip(labels, Xi):
        if np.any(row != 0):
            print(f"{lab:10s}" + "".join(f"{c:14.6g}" for c in row))


def sparse_galerkin(Xi, polyorder, usesine=False):
    """Right-hand side of the identified model (``sparseGalerkin``)."""
    tt = terms(polyorder, usesine)
    return lambda t, x: (library(np.atleast_2d(x), tt) @ Xi)[0]


def simulate(rhs, x0, t, rtol=1e-10, atol=1e-10):
    sol = solve_ivp(rhs, (t[0], t[-1]), x0, t_eval=t, rtol=rtol, atol=atol, method="DOP853")
    return sol.y.T


def central_difference(x, dt):
    """Fourth-order central difference on interior points (drops two samples at each end)."""
    return (-x[4:] + 8 * x[3:-1] - 8 * x[1:-3] + x[:-4]) / (12 * dt)


def color_line3(ax, x, y, z, c, lw=1.5, cmap="jet"):
    """A 3-D line coloured by ``c`` (``color_line3``)."""
    pts = np.column_stack([x, y, z])
    seg = np.stack([pts[:-1], pts[1:]], axis=1)
    lc = Line3DCollection(seg, cmap=cmap, linewidths=lw)
    lc.set_array(np.asarray(c)[:-1])
    ax.add_collection(lc)
    ax.set(xlim=(np.min(x), np.max(x)), ylim=(np.min(y), np.max(y)), zlim=(np.min(z), np.max(z)))
    return lc


def color_line(ax, x, y, c, lw=1.5, cmap="jet"):
    pts = np.column_stack([x, y])
    lc = LineCollection(np.stack([pts[:-1], pts[1:]], axis=1), cmap=cmap, linewidths=lw)
    lc.set_array(np.asarray(c)[:-1])
    ax.add_collection(lc)
    ax.autoscale()
    return lc


__all__ = ["stls", "terms", "pool_data", "print_model", "sparse_galerkin", "simulate", "central_difference",
           "color_line3", "color_line"]
