import numpy as np
import pytest

from unicodes import controls as c

A4 = np.array([[-0.0284, 6.6202, 0, -32.1890], [-0.0002, -0.9542, 0.9942, -0.012],
               [0.0013, -7.3610, -1.3895, 0.0005], [0, 0, 1, 0]])
B4 = np.array([[0], [-0.0748], [-15.0071], [0]])


def test_modes_short_period_and_phugoid():
    m = c.modes(A4)
    assert m[0].damping_ratio == pytest.approx(m[1].damping_ratio)
    assert m[0].natural_frequency > 1 > m[-1].natural_frequency  # short period fast, phugoid slow
    assert all(mode.eigenvalue.real < 0 for mode in m)


def test_transfer_function_dc_gain_matches_state_space():
    num, den = c.transfer_functions(A4, B4)
    dc_tf = num[:, -1] / den[-1]
    dc_ss = -np.linalg.solve(A4, B4).ravel()
    np.testing.assert_allclose(dc_tf, dc_ss, rtol=1e-6, atol=1e-10)


def test_doublet_shape():
    t = np.arange(0, 5, 0.01)
    u = c.doublet(t, 2, 0.5)
    assert u[np.searchsorted(t, 2.2)] == 1 and u[np.searchsorted(t, 2.7)] == -1 and u[0] == 0
    assert c.wrap_to_pi(3 * np.pi) == pytest.approx(-np.pi)


def test_mpc_tracks_pitch_within_bounds():
    C = np.array([[0, 0, 0, 1.0]])
    mpc = c.LinearMPC(A4, B4, C, dt=0.05, prediction_horizon=40, control_horizon=5, w_y=1.0, w_du=0.01,
                      u_min=-np.deg2rad(20), u_max=np.deg2rad(20))
    X, U, Y = mpc.simulate(np.zeros(4), np.deg2rad(5), 200)
    assert Y[-1, 0] == pytest.approx(np.deg2rad(5), abs=np.deg2rad(0.3))
    assert np.all(np.abs(U) <= np.deg2rad(20) + 1e-9)
