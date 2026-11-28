"""Particle swarm optimisation (AE 725 ``particle_swarm_optimization``)."""

from __future__ import annotations

import numpy as np


def particle_swarm_minimize(f, lower, upper, n_particles=50, max_iterations=100, inertia=0.5, c1=1.5, c2=1.5,
                            nonnegative=False, vectorized=False, rng=None):
    """Global-best PSO. Returns ``(best_x, best_f)``.

    ``f`` takes one point, or with ``vectorized=True`` an ``(n_dims,
    n_particles)`` array. Fixes from the MATLAB: personal bests were updated
    with a logical mask over the whole matrix (only the first row changed),
    were re-evaluated every iteration instead of stored, and the initial
    global best took a row instead of a column.
    """
    rng = np.random.default_rng(rng)
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    n = lower.size

    def evaluate(p):
        if vectorized:
            return np.asarray(f(p), float)
        return np.array([f(p[:, i]) for i in range(p.shape[1])], dtype=float)

    pos = lower[:, None] + rng.random((n, n_particles)) * (upper - lower)[:, None]
    vel = rng.random((n, n_particles))
    pbest, pbest_f = pos.copy(), evaluate(pos)
    i = int(np.argmin(pbest_f))
    gbest, gbest_f = pbest[:, i].copy(), pbest_f[i]
    for _ in range(max_iterations):
        vel = (inertia * vel + c1 * rng.random((n, n_particles)) * (pbest - pos)
               + c2 * rng.random((n, n_particles)) * (gbest[:, None] - pos))
        pos = pos + vel
        if nonnegative:
            pos = np.maximum(pos, 0) + 1e-12
        fit = evaluate(pos)
        better = fit < pbest_f
        pbest[:, better], pbest_f[better] = pos[:, better], fit[better]
        i = int(np.argmin(pbest_f))
        if pbest_f[i] < gbest_f:
            gbest, gbest_f = pbest[:, i].copy(), pbest_f[i]
    return gbest, float(gbest_f)
