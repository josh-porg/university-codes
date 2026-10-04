import numpy as np

from unicodes.optimize import (
    constrained_steepest_descent,
    golden_section_step,
    inexact_step_size,
    numeric_constraints,
    particle_swarm_minimize,
    qp_direction,
)
from unicodes.structures import WingBox, assemble_beam, displacement_sensitivity


def test_line_searches():
    assert abs(golden_section_step(lambda x: (x - 3) ** 2, 0.0, 1.0, tol=1e-8) - 3) < 1e-6
    t = inexact_step_size(lambda x: x @ x, np.array([2.0, 0]), np.array([-2.0, 0]))
    assert t == 1.0


def test_qp_direction_matches_kkt():
    # min c.d + d.d/2 s.t. d1 + d2 <= -1 -> d = -(c + u [1, 1])
    d, u = qp_direction([1.0, 0.0], [[1.0, 1.0]], [-3.0])
    assert np.isclose(d.sum(), -3, atol=1e-8)
    assert np.allclose(d, -(np.array([1.0, 0]) + u[0]), atol=1e-8)


def test_csd_beam_problem():
    def fun(x):
        return x[0] * x[1], np.array([x[1], x[0]])

    def g(x):
        b, d = x
        return np.array([2.4e8 / (b * d**2) - 10, 4.5e5 / (2 * b * d) - 2, d - 2 * b])

    r = constrained_steepest_descent(fun, numeric_constraints(g, 1e-6), [250, 300], max_iterations=2000)
    assert np.isclose(r.fun, 112500, rtol=1e-3)
    assert g(r.x).max() < 1e-3


def test_pso_sphere():
    x, f = particle_swarm_minimize(lambda x: np.sum((x - 1) ** 2), -5 * np.ones(3), 5 * np.ones(3), n_particles=40,
                                   max_iterations=200, rng=0)
    assert f < 1e-6 and np.allclose(x, 1, atol=1e-3)


def test_wingbox_equilibrium():
    wb = WingBox()
    x = np.array([1.2, 0.8, 0.8, 1.2, 0.05, 0.05, 0.05, 0.05])
    B = wb.boom_areas(x)
    y, z = wb.boom_coordinates()
    s = wb.stresses(x)
    # axial forces balance and reproduce the bending moment
    sigma = np.zeros(10)
    sigma[[0, 5, 6, 9]] = s.cap_axial
    sigma[[1, 2, 3, 4, 7, 8]] = s.stringer_axial
    assert abs(B @ sigma) < 1e-6 * wb.M_z
    assert np.isclose(-(B * sigma) @ (y - B @ y / B.sum()), wb.M_z)
    assert wb.constraints(x, buckling=True).size == 13 + 8
    assert wb.constraints(np.r_[x, [0.08] * 6]).size == 14 + 14


def test_beam_sensitivity_direct_equals_adjoint_and_fd():
    K = assemble_beam(4.0, [2.0, 1.0], 1.0)
    F = np.zeros(6)
    F[4] = 36
    free = np.arange(2, 6)
    tip = np.array([0, 0, 1, 0])
    dK = assemble_beam(4.0, [1.0, 0.0], 1.0)
    direct = displacement_sensitivity(K, dK, F, free, tip)
    adjoint = displacement_sensitivity(K, dK, F, free, tip, adjoint=True)
    h = 1e-6
    w = lambda I1: np.linalg.solve(assemble_beam(4.0, [I1, 1.0], 1.0)[2:, 2:], F[2:])[2]  # noqa: E731
    assert np.isclose(direct, adjoint)
    assert np.isclose(direct, (w(2 + h) - w(2 - h)) / (2 * h), rtol=1e-5)
