"""Real-coded genetic algorithm (AE 725 ``genetic_algorithm_optimization``/``GeneticAlgorithm``,
and the GA loop of the AE 722 ``AETHER_Ben_6DOF_with_optimizer*`` scripts).

Selection is roulette-wheel on ``1 / (1 + f - f_min)``, crossover is one-point
on a fraction of the mating pool, mutation either adds N(0, sigma) noise to
genes (AE 725) or redraws a gene uniformly inside the bounds (AE 722), and
the best ``n_elites`` survive unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class GAResult:
    x: np.ndarray
    fun: float
    history: list = field(default_factory=list)  # best value per generation
    population: np.ndarray | None = None


def genetic_minimize(
    f,
    lower,
    upper,
    population_size=50,
    generations=100,
    crossover_rate=0.8,
    mutation_rate=0.1,
    n_elites=1,
    mutation="gaussian",
    sigma=1.0,
    clip_to_bounds=False,
    positive=False,
    vectorized=False,
    rng=None,
    callback=None,
) -> GAResult:
    """Minimise ``f`` with a genetic algorithm.

    ``lower``/``upper`` bound the initial population (and mutation when
    ``mutation="uniform"``; with ``clip_to_bounds`` every generation is
    clipped). ``positive=True`` reproduces the AE 725 wing-box version that
    forces genes to stay above zero. With ``vectorized=True`` ``f`` gets the
    whole population as a ``(n_genes, population_size)`` array.
    """
    rng = np.random.default_rng(rng)
    lower, upper = np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)
    n = lower.size
    pop = lower[:, None] + rng.random((n, population_size)) * (upper - lower)[:, None]

    def evaluate(p):
        if vectorized:
            return np.asarray(f(p), dtype=float)
        return np.array([f(p[:, i]) for i in range(p.shape[1])], dtype=float)

    history = []
    fit = evaluate(pop)
    for gen in range(generations):
        order = np.argsort(fit)
        history.append(float(fit[order[0]]))
        if callback:
            callback(gen, pop[:, order[0]], fit[order[0]])
        w = 1 / (1 + fit - fit.min())
        pool = pop[:, rng.choice(population_size, population_size, p=w / w.sum())]

        new = pool.copy()
        n_cross = int(round(crossover_rate * population_size / 2)) * 2
        if n > 1 and n_cross:
            idx = rng.permutation(population_size)[:n_cross]
            for a, b in zip(idx[::2], idx[1::2]):
                cp = rng.integers(1, n)
                new[cp:, a], new[cp:, b] = pool[cp:, b].copy(), pool[cp:, a].copy()

        mask = rng.random(new.shape) < mutation_rate
        if mutation == "gaussian":
            new = new + mask * rng.normal(0, sigma, new.shape)
        else:
            redraw = lower[:, None] + rng.random(new.shape) * (upper - lower)[:, None]
            new = np.where(mask, redraw, new)
        if positive:
            new = np.maximum(new, 0) + 1e-12
        if clip_to_bounds:
            new = np.clip(new, lower[:, None], upper[:, None])

        new[:, :n_elites] = pop[:, order[:n_elites]]
        pop = new
        fit = evaluate(pop)

    best = int(np.argmin(fit))
    history.append(float(fit[best]))
    return GAResult(pop[:, best].copy(), float(fit[best]), history, pop)
