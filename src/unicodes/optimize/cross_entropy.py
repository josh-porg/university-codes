"""Cross-entropy method, as used for the propeller and guidance optimisations."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np


@dataclass
class CrossEntropyResult:
    x: np.ndarray
    fun: float
    history: list[float]
    n_iterations: int


def cross_entropy_maximize(
    objective: Callable[[np.ndarray], float],
    x0,
    rel_std,
    lower=None,
    upper=None,
    n_samples: int = 300,
    elite_fraction: float = 0.6,
    max_iterations: int = 18,
    rtol: float = 2e-7,
    rng=None,
    map_fn=map,
) -> CrossEntropyResult:
    """Maximise ``objective`` by iteratively resampling around the best candidates.

    Each iteration draws ``x = mean * (1 + rel_std * N(0, 1))`` (clipped
    to the bounds), keeps the best ``elite_fraction``, moves the mean to the
    best sample and sets ``rel_std`` from the spread of the elite. Stops
    when the best value improves by less than ``rtol`` (relative).

    ``map_fn`` lets you evaluate samples in parallel, e.g.
    ``concurrent.futures.ProcessPoolExecutor().map``.
    """
    rng = np.random.default_rng(rng)
    mean = np.asarray(x0, dtype=float)
    rel_std = np.broadcast_to(np.asarray(rel_std, dtype=float), mean.shape).copy()
    lower = np.full_like(mean, -np.inf) if lower is None else np.asarray(lower, dtype=float)
    upper = np.full_like(mean, np.inf) if upper is None else np.asarray(upper, dtype=float)

    best_x, best_f = mean.copy(), float(objective(mean))
    history = [best_f]
    it = 0
    for it in range(1, max_iterations + 1):
        samples = mean * (1 + rel_std * rng.standard_normal((n_samples, mean.size)))
        samples = np.clip(samples, lower, upper)
        values = np.array(list(map_fn(objective, samples)), dtype=float)

        order = np.argsort(values)[::-1]
        elite = samples[order[: max(2, round(elite_fraction * n_samples))]]
        previous = mean
        mean = samples[order[0]]
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_std = np.nan_to_num(np.abs(elite.std(axis=0, ddof=1) / previous))

        last = best_f
        if values[order[0]] > best_f:
            best_x, best_f = mean.copy(), float(values[order[0]])
        history.append(best_f)
        if last > 0 and (best_f - last) / last < rtol:
            break
    return CrossEntropyResult(best_x, best_f, history, it)
