import numpy as np
import pytest

from unicodes import orbital

R0 = np.array([1131.34, -2282.343, 6672.423]) * 1e3  # AE 360 HW9 initial state
V0 = np.array([-5.64305, 4.30333, 2.42879]) * 1e3


def test_coe_roundtrip():
    coe = orbital.rv_to_coe(R0, V0)
    r, v = orbital.coe_to_rv(coe)
    np.testing.assert_allclose(r, R0, rtol=1e-9)
    np.testing.assert_allclose(v, V0, rtol=1e-9)


def test_kepler_solution_satisfies_equation():
    for M, e in [(0.3, 0.1), (-2.0, 0.7), (3.1, 0.95)]:
        E = orbital.solve_kepler(M, e)
        assert E - e * np.sin(E) == pytest.approx(M, abs=1e-12)


def test_analytic_and_numerical_propagation_agree():
    dt = 24 * 3600.0
    r_k, v_k = orbital.propagate_kepler(R0, V0, dt)
    sol = orbital.propagate(R0, V0, (0, dt))
    np.testing.assert_allclose(sol.y[:3, -1], r_k, atol=1.0)  # metres after one day
    E0 = orbital.specific_energy(R0, V0)
    assert orbital.specific_energy(r_k, v_k) == pytest.approx(E0, rel=1e-10)


def test_euler_integration_drifts_in_energy():
    _, r, v = orbital.propagate_euler(R0, V0, duration=3600, dt=1.0)
    E = orbital.specific_energy(r, v)
    assert abs(E[-1] - E[0]) / abs(E[0]) > 1e-6  # the homework's point: Euler is not conservative
