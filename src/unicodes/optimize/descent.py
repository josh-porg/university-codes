"""Line searches and constrained steepest descent (AE 725, Math 796).

Ports of ``golden_search``, ``inexactStepSize``, ``DescentFunction``,
``constrained_steepest_descent``, ``continuous_simulated_anealing`` and the
Metropolis loop of ``Wingbox_SA_V_0``.

The constrained methods follow Arora's CSD algorithm: at each iterate the
quadratic subproblem ``min c.d + 0.5 d.d  s.t.  A d <= -g`` gives the search
direction ``d`` and the multipliers ``u``; the penalty parameter
``R = max(R, sum u)`` defines the descent function ``phi = f + R V`` with
``V = max(0, g_i)``, and the step length comes from an Armijo search on
``phi``. ``fun`` returns ``(f, grad_f)`` and ``constraints`` returns
``(g, dg/dx)`` with ``g <= 0`` feasible, as in the MATLAB.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from ..numerics import jacobian

GOLDEN = (1 + np.sqrt(5)) / 2


def golden_section_step(f, x, d, tol=1e-6, delta=0.1):
    """Step length minimising ``f(x + alpha d)`` by golden-section search (``golden_search``).

    Phase one brackets the minimum by growing the step by ``delta * phi^q``;
    phase two shrinks the bracket to ``tol``. Fixes from the MATLAB: the
    delta-reduction loop evaluated ``f(x + delta)`` instead of
    ``f(x + delta d)``, and the tie branch moved ``alpha_b`` the wrong way.
    """
    x, d = np.asarray(x, float), np.asarray(d, float)

    def phi(a):
        return f(x + a * d)

    f0 = phi(0.0)
    while phi(delta) > f0 and delta > 1e-12:
        delta *= 0.1
    q = 0
    alpha = delta
    f_old, f_new = f0, phi(alpha)
    while f_new < f_old:
        f_old = f_new
        q += 1
        alpha += delta * GOLDEN**q
        f_new = phi(alpha)
    upper = alpha
    lower = max(alpha - delta * GOLDEN**q - delta * GOLDEN ** (q - 1), 0.0)

    tau = GOLDEN - 1
    a = upper - tau * (upper - lower)
    b = lower + tau * (upper - lower)
    fa, fb = phi(a), phi(b)
    while upper - lower > tol:
        if fa < fb:
            upper, b, fb = b, a, fa
            a = upper - tau * (upper - lower)
            fa = phi(a)
        elif fa > fb:
            lower, a, fa = a, b, fb
            b = lower + tau * (upper - lower)
            fb = phi(b)
        else:
            lower, upper = a, b
            a = upper - tau * (upper - lower)
            b = lower + tau * (upper - lower)
            fa, fb = phi(a), phi(b)
    return (lower + upper) / 2


def inexact_step_size(phi, x, d, gamma=0.5, mu=0.5, max_trials=60):
    """Armijo backtracking (``inexactStepSize``): the first ``t = mu^j`` with
    ``phi(x + t d) <= phi(x) - t gamma |d|^2``."""
    x, d = np.asarray(x, float), np.asarray(d, float)
    phi0 = phi(x)
    beta = gamma * float(np.sum(d * d))
    t = 1.0
    for _ in range(max_trials):
        if phi(x + t * d) <= phi0 - t * beta:
            break
        t *= mu
    return t


def descent_function(fun, constraints, R_old, multipliers):
    """Penalty descent function ``phi(x) = f(x) + R max(0, g(x))`` and the updated ``R`` (``DescentFunction``)."""
    R = max(R_old, float(np.sum(multipliers)))

    def phi(x):
        return fun(x)[0] + R * max(0.0, float(np.max(constraints(x)[0])))

    return phi, R


def qp_direction(c, A, b):
    """Solve ``min c.d + 0.5 d.d  s.t.  A d <= b`` (the CSD subproblem, ``quadprog`` with H = I).

    Solved through its dual ``min_{u >= 0} 0.5 |c + A^T u|^2 + b.u``, which
    gives ``d = -(c + A^T u)`` and the multipliers ``u`` directly.
    """
    c, A, b = np.asarray(c, float), np.atleast_2d(np.asarray(A, float)), np.asarray(b, float)
    if A.size == 0:
        return -c, np.zeros(0)

    def dual(u):
        r = c + A.T @ u
        return 0.5 * r @ r + b @ u, A @ r + b

    res = minimize(dual, np.zeros(A.shape[0]), jac=True, method="L-BFGS-B",
                   bounds=[(0, None)] * A.shape[0], options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 10000})
    u = res.x
    return -(c + A.T @ u), u


def numeric_constraints(g, step=1e-3):
    """Wrap ``g(x)`` to return ``(g, dg/dx)`` with a relative central-difference Jacobian
    (the drivers' ``inequalConstraintFunction`` + ``computeJacobian``)."""

    def wrapped(x):
        return np.atleast_1d(g(x)), jacobian(g, x, step)

    return wrapped


@dataclass
class CSDResult:
    x: np.ndarray
    fun: float
    iterations: int
    converged: bool
    x_history: list = field(default_factory=list)
    f_history: list = field(default_factory=list)
    d_history: list = field(default_factory=list)
    alpha_history: list = field(default_factory=list)
    T_history: list = field(default_factory=list)


def constrained_steepest_descent(fun, constraints, x0, max_iterations=100, eps1=1e-3, eps2=1e-3,
                                 R=1.0, gamma=0.5, mu=0.5, T0=0.0, freeze_time=-1, cooling="log", beta=0.99,
                                 rng=None):
    """Constrained steepest descent, optionally with continuous simulated annealing.

    ``T0 = 0`` (or ``freeze_time = -1``) is plain CSD
    (``constrained_steepest_descent``). Otherwise each step adds
    ``T U(-1, 1)`` noise until iteration ``freeze_time``
    (``continuous_simulated_anealing``). ``cooling`` is ``"log"``
    (``T0 / log(2 + k)``), ``"linear"`` (``T0 - beta k``) or
    ``"exponential"`` (``T0 beta^k``), the schedules tried in
    ``GlobalOptimization_*_Test``. Stops when the maximum constraint is below
    ``eps1`` and ``|d| < eps2``.
    """
    schedules = {
        "log": lambda k: T0 / np.log(2 + k),
        "linear": lambda k: max(T0 - beta * k, 0.0),
        "exponential": lambda k: T0 * beta**k,
    }
    schedule = schedules[cooling]
    rng = np.random.default_rng(rng)
    x = np.asarray(x0, dtype=float).copy()
    out = CSDResult(x, np.inf, 0, False)
    for k in range(max_iterations):
        x = np.where(x == 0, 1e-6, x)  # the SA version nudged zero variables off the boundary
        f, c = fun(x)
        g, A = constraints(x)
        g = np.atleast_1d(g)
        out.x_history.append(x.copy())
        out.f_history.append(float(f))
        d, u = qp_direction(c, A, -g)
        if np.max(g) <= eps1 and np.linalg.norm(d) <= eps2:
            out.converged = True
            break
        phi, R = descent_function(fun, constraints, R, u)
        alpha = inexact_step_size(phi, x, d, gamma, mu)
        T = schedule(k) if k <= freeze_time else 0.0
        x = x + alpha * d + T * rng.uniform(-1, 1, x.shape)
        out.d_history.append(d)
        out.alpha_history.append(alpha)
        out.T_history.append(T)
        out.iterations = k + 1
    out.x = x
    out.fun = float(fun(x)[0])
    return out


def simulated_annealing_minimize(f, x0, T0=100.0, max_iterations=100_000, T_min=1e-6, step=1.0, bias=0.0,
                                 lower=None, rng=None):
    """Metropolis simulated annealing with logarithmic cooling ``T0 / log(2 + k)`` (``Wingbox_SA_V_0``).

    Trial points are ``x + step N(0, 1) + bias`` (the MATLAB used ``bias=1e-2``)
    clipped at ``lower``. Returns ``(best_x, best_f)``.
    """
    rng = np.random.default_rng(rng)
    x = np.asarray(x0, float).copy()
    best_x, best_f = x.copy(), f(x)
    for k in range(1, max_iterations + 1):
        x_new = x + step * rng.standard_normal(x.shape) + bias
        if lower is not None:
            x_new = np.maximum(x_new, lower)
        f_new = f(x_new)
        # The MATLAB compared to the best value instead of the current one; kept.
        df = f_new - best_f
        T = T0 / np.log(2 + k)
        if df < 0 or rng.random() < np.exp(-df / T):
            x = x_new
            if f_new < best_f:
                best_x, best_f = x_new.copy(), f_new
        if T < T_min:
            break
    return best_x, best_f
