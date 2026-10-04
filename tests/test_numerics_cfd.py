import numpy as np
import pytest

from unicodes.cfd.advection import SCHEMES, advect
from unicodes.cfd.flux_reconstruction import fr_operators
from unicodes.cfd.mesh import resample_by_arclength, ruled_mesh
from unicodes.cfd.nozzle import NozzleFlow, characteristic_state, nozzle_area, solve_nozzle
from unicodes.ode import INTEGRATORS


@pytest.mark.parametrize("name, order", [("explicit_euler", 1), ("modified_euler", 2), ("ssp_rk3", 3), ("rk4", 4)])
def test_integrator_order(name, order):
    f = INTEGRATORS[name]
    errs = [abs(f(lambda t, y: -y, (0, 1), dt, 1.0)[1][-1] - np.exp(-1)) for dt in (0.02, 0.01)]
    assert abs(np.log2(errs[0] / errs[1]) - order) < 0.15


def test_integrators_handle_vectors():
    t, y = INTEGRATORS["rk4"](lambda t, y: np.array([y[1], -y[0]]), (0, np.pi), np.pi / 300, [0.0, 1.0])
    assert np.allclose(y[-1], [0, -1], atol=1e-6)


def test_advection_conserves_mass():
    x = np.linspace(0, 10, 80, endpoint=False)
    u0 = np.where(x < 2 * np.pi, 1 - np.cos(x), 0)
    for name in ("upwind1", "upwind2_ssprk2"):
        assert np.isclose(advect(u0, 0.25, 100, name).sum(), u0.sum())
    assert set(SCHEMES) == {"upwind1", "central2_euler", "upwind2_ssprk2"}


@pytest.mark.parametrize("case", [1, 2])
def test_exact_nozzle_mass_flow(case):
    flow = NozzleFlow(case)
    x = np.linspace(-4, 4, 41)
    Q = flow(x)
    mdot = Q[:, 1] * nozzle_area(x)
    assert np.allclose(mdot, mdot[0], rtol=1e-6)
    if case == 2:
        assert 0 < flow.x_shock < 4


def test_characteristic_state_is_identity_for_equal_states():
    Q = NozzleFlow(1)(0.5)
    assert np.allclose(characteristic_state(Q, Q), Q)


def test_nozzle_solver_subsonic():
    flow = NozzleFlow(1)
    r = solve_nozzle(flow, 60, order=2, limiter="vanleer", tol=1e-5)
    assert r.converged
    assert np.sqrt(np.mean((r.Q[:, 0] - flow(r.x)[:, 0]) ** 2)) < 0.05


def test_fr_operators():
    for k in range(5):
        op = fr_operators(k)
        assert np.allclose(op.D @ op.nodes**k, k * op.nodes ** max(k - 1, 0) if k else 0)
        assert np.allclose(op.interp @ op.nodes**k, [(-1) ** k, 1])
        assert np.allclose(op.dg_left.sum(), -op.dg_right.sum())


def test_ruled_mesh():
    inner = np.column_stack([np.linspace(0, 1, 5), np.zeros(5)])
    outer = np.column_stack([np.linspace(0, 2, 3), np.ones(3)])
    nodes = ruled_mesh(inner, outer, 11, 4)
    assert nodes.shape == (11, 4, 2)
    assert np.allclose(nodes[:, -1, 0], np.linspace(0, 2, 11))
    assert np.allclose(resample_by_arclength([[0, 0], [1, 0], [1, 1]], 3), [[0, 0], [1, 0], [1, 1]])
