"""Optimisers used across projects (AE 722, AE 725, Math 796)."""

from .cross_entropy import cross_entropy_maximize
from .descent import (
    CSDResult,
    constrained_steepest_descent,
    descent_function,
    golden_section_step,
    inexact_step_size,
    numeric_constraints,
    qp_direction,
    simulated_annealing_minimize,
)
from .genetic import GAResult, genetic_minimize
from .swarm import particle_swarm_minimize

__all__ = [
    "CSDResult",
    "GAResult",
    "constrained_steepest_descent",
    "cross_entropy_maximize",
    "descent_function",
    "genetic_minimize",
    "golden_section_step",
    "inexact_step_size",
    "numeric_constraints",
    "particle_swarm_minimize",
    "qp_direction",
    "simulated_annealing_minimize",
]
