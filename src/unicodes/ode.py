"""Fixed-step explicit ODE integrators (AE 746 ``explicitEuler``, ``modifiedEuler``, ``SSP_RK3``, ``RK4``;
the AE 722 PAH ``RK4``).

All take ``f(t, y)``, ``t_span = (t0, t1)``, step ``dt`` and ``y0`` (scalar or
array) and return ``(t, y)`` with ``y[k]`` the state at ``t[k]``. The MATLAB
versions stored only the first state component; vectors work here.
"""

from __future__ import annotations

import numpy as np


def _grid(t_span, dt, y0):
    n = int(round((t_span[1] - t_span[0]) / dt))
    t = t_span[0] + dt * np.arange(n + 1)
    y0 = np.asarray(y0, dtype=float)
    y = np.empty((n + 1,) + y0.shape)
    y[0] = y0
    return t, y


def explicit_euler(f, t_span, dt, y0):
    t, y = _grid(t_span, dt, y0)
    for i in range(len(t) - 1):
        y[i + 1] = y[i] + dt * np.asarray(f(t[i], y[i]))
    return t, y


def modified_euler(f, t_span, dt, y0):
    """Heun's method: Euler predictor, trapezoidal corrector."""
    t, y = _grid(t_span, dt, y0)
    for i in range(len(t) - 1):
        k1 = np.asarray(f(t[i], y[i]))
        k2 = np.asarray(f(t[i + 1], y[i] + dt * k1))
        y[i + 1] = y[i] + dt / 2 * (k1 + k2)
    return t, y


def ssp_rk3(f, t_span, dt, y0):
    """Shu-Osher strong-stability-preserving third-order Runge-Kutta."""
    t, y = _grid(t_span, dt, y0)
    for i in range(len(t) - 1):
        y1 = y[i] + dt * np.asarray(f(t[i], y[i]))
        y2 = 0.75 * y[i] + 0.25 * (y1 + dt * np.asarray(f(t[i] + dt, y1)))
        y[i + 1] = y[i] / 3 + 2 / 3 * (y2 + dt * np.asarray(f(t[i] + dt / 2, y2)))
    return t, y


def rk4(f, t_span, dt, y0):
    """Classical fourth-order Runge-Kutta."""
    t, y = _grid(t_span, dt, y0)
    for i in range(len(t) - 1):
        k1 = dt * np.asarray(f(t[i], y[i]))
        k2 = dt * np.asarray(f(t[i] + dt / 2, y[i] + k1 / 2))
        k3 = dt * np.asarray(f(t[i] + dt / 2, y[i] + k2 / 2))
        k4 = dt * np.asarray(f(t[i] + dt, y[i] + k3))
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return t, y


INTEGRATORS = {"explicit_euler": explicit_euler, "modified_euler": modified_euler, "ssp_rk3": ssp_rk3, "rk4": rk4}
