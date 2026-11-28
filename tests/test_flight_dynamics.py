import numpy as np
import pytest

from unicodes.flight_dynamics import ride_quality as rq
from unicodes.flight_dynamics import sixdof as s6
from unicodes.flight_dynamics import turbulence as tb
from unicodes.numerics import jacobian
from unicodes.optimize import genetic_minimize


def test_sixdof_level_flight_is_steady_without_inputs():
    ac = s6.cessna_182()
    x0 = s6.initial_state(ac, theta=0.0)
    xdot = s6.derivatives(0, x0, np.zeros(4), ac)
    # trimmed at alpha = 0 with C_L1 lift and constant thrust: only small residual accelerations
    assert abs(xdot[2]) < 1e-12 and abs(xdot[3]) < 1e-12  # no lateral motion
    assert xdot[9] == pytest.approx(ac.V_trim)


def test_sixdof_modes_are_stable_and_include_phugoid():
    ac = s6.cessna_182()
    J = s6.linearize(ac, s6.initial_state(ac, theta=0.0), np.zeros(4))
    long = [0, 1, 4, 7]  # V, alpha, theta, Q
    lam = np.linalg.eigvals(J[np.ix_(long, long)])
    assert np.all(lam.real < 0)
    assert np.min(np.abs(lam.imag[lam.imag > 0])) < 0.5  # phugoid: slow oscillation


def test_flap_spring_returns_to_neutral():
    f = s6.PAHFlap(0.25, [0, 0, 0], [0, 1, 0], [-0.25, 0, 0], [0, 0.5, 0], I_h=1.0, k=5.0, c_r=1.0,
                   delta_0=0.1, C_h_alpha=0.0, C_h_delta=0.0)
    assert f.hinge_acceleration(0.1, 0.0, 1.2, np.array([30, 0, 0]), (0, 0, 0)) == pytest.approx(0)
    assert f.hinge_acceleration(0.0, 0.0, 1.2, np.array([30, 0, 0]), (0, 0, 0)) == pytest.approx(0.5)
    f.c_c = 1.0  # Coulomb friction larger than the spring moment holds it
    assert f.hinge_acceleration(0.0, 0.0, 1.2, np.array([30, 0, 0]), (0, 0, 0)) == 0.0


def test_dryden_rms_matches_sigma():
    ug, vg, wg = tb.dryden_gusts(200_000, 0.1, 200.0, 5000.0, rng=0)
    sigma = 15 * 0.3048
    for g in (ug, vg, wg):
        assert np.std(g) == pytest.approx(sigma, rel=0.05)


def test_iso_weighting_shape_and_comfort():
    num, den = rq.iso2631_weighting("Wk")
    H = lambda f: abs(np.polyval(num, 2j * np.pi * f) / np.polyval(den, 2j * np.pi * f))  # noqa: E731
    assert H(6.3) == pytest.approx(1.05, abs=0.05)  # ISO 2631-1 Table 3: Wk ~ 1.05 near 5-8 Hz
    assert H(0.1) < 0.1 and H(400) < 0.1
    assert rq.comfort_level(0.2, 0.3) == "Not uncomfortable"
    assert rq.comfort_level(0.2, 1.2) == "Fairly uncomfortable"
    t = np.arange(0, 10, 0.01)
    r = rq.ride_quality(t, 0 * t, 0 * t, np.sin(2 * np.pi * 5 * t), 100)
    assert r.rms[2] == pytest.approx(1.05 / np.sqrt(2), rel=0.1)


def test_jacobian_and_ga():
    J = jacobian(lambda x: np.array([x[0] ** 2, x[0] * x[1]]), [3.0, 2.0], relative=False)
    np.testing.assert_allclose(J, [[6, 0], [2, 3]], atol=1e-6)
    res = genetic_minimize(lambda x: np.sum((x - 1) ** 2), [-5, -5], [5, 5], population_size=60,
                           generations=60, mutation_rate=0.2, sigma=0.3, n_elites=2, rng=0)
    assert res.fun < 0.05
    assert all(np.diff(res.history) <= 1e-12)  # elitism: best never gets worse


# --- regression tests: gust enters the aerodynamics, not the kinematics; flap roll arm is y/b ---


def _c182_flaps(ac, k=3.5, zeta=0.7):
    c_f = 0.4 * ac.c_bar
    q = 0.5 * ac.rho * ac.V_trim**2
    I_h, C_ha, C_hd = 0.0674, -0.225, -0.370
    wn = np.sqrt((k + q * c_f**2 * abs(C_hd)) / I_h)
    return s6.wing_flaps(ac.b / 2, c_f, 0.748, I_h=I_h, k=k, c_r=2 * zeta * I_h * wn, delta_0=0.0,
                         C_h_alpha=C_ha, C_h_delta=C_hd,
                         delta_min=-np.deg2rad(30), delta_max=np.deg2rad(30))


def test_updraft_changes_rigid_lift_but_not_kinematics():
    ac = s6.cessna_182()
    x0 = s6.initial_state(ac, theta=0.0)
    calm = s6.derivatives(0, x0, np.zeros(4), ac)
    w_up = -10.0  # ft/s; w_g < 0 is an updraft
    gusty = s6.derivatives(0, x0, np.zeros(4), ac, gust=np.array([0.0, 0.0, w_up]))
    # inertial position rates do not move with the air mass
    np.testing.assert_allclose(gusty[9:12], calm[9:12], atol=1e-12)
    # the wing feels the air-relative angle of attack: extra upward (-z) specific force
    dfz = s6.specific_force(x0, np.zeros(4), ac, gust=np.array([0.0, 0.0, w_up]))[2] - s6.specific_force(x0, np.zeros(4), ac)[2]
    d_alpha = np.arctan2(-w_up, ac.V_trim)
    q_a = 0.5 * ac.rho * (ac.V_trim**2 + w_up**2)
    expected = -q_a * ac.S * ac.derivatives.C_L_alpha * d_alpha / ac.mass
    assert dfz < 0
    assert abs(dfz / expected - 1) < 0.1
    # and the gust drives the short-period response (alpha_dot, Q_dot change)
    assert abs(gusty[1] - calm[1]) > 1e-3 and abs(gusty[7] - calm[7]) > 1e-3


def test_flap_rolling_moment_uses_nondimensional_arm():
    ac = s6.cessna_182()
    flaps = _c182_flaps(ac)
    x = s6.initial_state(ac, theta=0.0, flaps=flaps)
    x[16], x[18] = 0.05, -0.05  # anti-symmetric flap deflection
    Pdot = s6.derivatives(0, x, np.zeros(4), ac, flaps)[6]
    q = 0.5 * ac.rho * ac.V_trim**2
    d = ac.derivatives
    C_l = -sum(0.5 * d.C_L_alpha * f.tau * x[16 + 2 * i] * f.r_cg[1] / ac.b for i, f in enumerate(flaps))
    np.testing.assert_allclose(Pdot, C_l * q * ac.S * ac.b / ac.I_xx, rtol=1e-6)


def test_two_flap_cessna_stays_bounded_in_turbulence():
    from unicodes.flight_dynamics.turbulence import dryden_gusts

    FT = 0.3048
    ac = s6.cessna_182()
    flaps = _c182_flaps(ac)
    dt, T = 0.02, 10.0
    t = np.arange(0, T + dt / 2, dt)
    ug, vg, wg = dryden_gusts(len(t), dt, ac.V_trim * FT, 5000 * FT, "severe", rng=0, sigma_high=15 * FT)
    gusts = np.column_stack([ug, vg, wg]) / FT
    x0 = s6.initial_state(ac, theta=0.0, flaps=flaps)
    res = s6.simulate(ac, t, x0, np.zeros((len(t), 4)), flaps, gusts, rtol=1e-4, atol=1e-6)
    assert np.all(np.isfinite(res.x))
    assert np.max(np.abs(res["phi"])) < np.deg2rad(30)
    assert np.max(np.abs(res["delta_f1"])) < np.deg2rad(10)
    # the rigid aircraft responds to the gust through the air-relative alpha
    rigid = s6.simulate(ac, t, s6.initial_state(ac, theta=0.0), np.zeros((len(t), 4)), (), gusts, rtol=1e-4, atol=1e-6)
    assert np.std(rigid["alpha"]) > np.deg2rad(0.2)
    assert rigid.specific_force.shape == (len(t), 3)
    assert np.std(rigid.specific_force[:, 2]) > 1.0  # ft/s^2
