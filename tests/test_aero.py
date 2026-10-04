import numpy as np
import pytest

from unicodes.aero import Airfoil, PropellerBlade, soaring
from unicodes.optimize import cross_entropy_maximize


def make_airfoil():
    d = np.deg2rad
    return Airfoil("test", "unit test", 2 * np.pi, d(-4), 0.44, d(8), d(0), 1.3, d(14), 1.6, d(18), 1.1)


def test_airfoil_spline_hits_points_and_peaks_at_clmax():
    af = make_airfoil()
    alpha, c_l = af.key_points()
    np.testing.assert_allclose(af.c_l(alpha), c_l, atol=1e-12)
    assert af.c_l_alpha_at(af.alpha_c_l_max) == pytest.approx(0, abs=1e-9)
    grid = np.linspace(af.alpha_0, af.alpha_post_stall, 500)
    assert af.c_l(grid).max() == pytest.approx(af.c_l_max, abs=1e-3)
    assert af.c_l(np.deg2rad(25)) == 0.0


def test_airfoil_file(tmp_path):
    path = tmp_path / "af.txt"
    make_airfoil().write_airfoil_file(path)
    assert "alpha_c_l_max [deg]=14.0" in path.read_text()


POLAR = dict(C_D0=0.012, A=22.8, e=0.9, wing_loading=280.0, rho=1.0834)


def test_turn_sink_tends_to_straight_glide():
    straight = soaring.sink_rate(1.0, **POLAR)
    turning = soaring.sink_rate_in_turn(1.0, **POLAR, r=1e6)
    assert turning == pytest.approx(straight, rel=1e-6)
    assert np.isnan(soaring.sink_rate_in_turn(1.0, **POLAR, r=5.0))


def test_maccready_speed_to_fly_increases_with_climb_rate():
    cl_weak, v_weak = soaring.optimal_glide_cl(1.0, **POLAR, C_L_max=1.5)
    cl_strong, v_strong = soaring.optimal_glide_cl(3.0, **POLAR, C_L_max=1.5)
    assert cl_strong < cl_weak  # fly faster between strong thermals
    assert v_strong > v_weak


def test_best_climb_in_thermal_is_positive_for_strong_thermal():
    climb, r = soaring.best_climb_in_thermal("A2", **POLAR, C_L_max=1.4)
    assert 0 < climb < 3.5
    assert r > 0


def test_propeller_matches_requested_power():
    blade = PropellerBlade(np.deg2rad(45), np.deg2rad(15), np.deg2rad(-42), -0.7, 0.1, radius=1.1)
    V, rho = 75 * 0.514444, 0.002133 * 515.379
    perf = blade.match_power(19.7e3, V=V, rho=rho)
    assert perf.aero_power == pytest.approx(19.7e3, rel=1e-6)
    assert 0.5 < perf.efficiency < 1.0


def test_cross_entropy_finds_maximum():
    res = cross_entropy_maximize(
        lambda x: -np.sum((x - np.array([2.0, -1.0])) ** 2), x0=[1.0, -0.5], rel_std=0.5,
        n_samples=200, elite_fraction=0.2, max_iterations=40, rtol=0, rng=0,
    )
    np.testing.assert_allclose(res.x, [2.0, -1.0], atol=1e-2)


def test_ferry_mission_matches_ae521_script():
    from unicodes.aero import sizing

    V_cruise = sizing.mach_to_kts(0.75, 968.1)
    V_divert = sizing.mach_to_kts(0.75, 1036.8)
    mission = [
        ("engine start", 0.99), ("taxi", 0.99), ("take-off", 0.995), ("climb", 0.985),
        ("cruise", sizing.breguet_range_fraction(3000, V_cruise, 0.34, 26)),
        ("descent", 0.995), ("landing", 1.0), ("climb 2", 0.995),
        ("divert", sizing.breguet_range_fraction(100, V_divert, 0.34, 26)),
        ("loiter", sizing.breguet_endurance_fraction(0.75, 0.34, 26)),
        ("descent 2", 0.995), ("landing 2", 1.0), ("shutdown", 0.995),
    ]
    W_TO = 130000 - 75000
    est = sizing.estimate_weights(W_TO, mission, payload=0, crew=700)
    M_ff = np.prod([f for _, f in mission])
    assert est.fuel_used == pytest.approx((1 - M_ff) * W_TO)
    assert est.empty_tentative == pytest.approx(W_TO - est.fuel_used - 700 - 0.01 * est.fuel_used)


def test_payload_drop_reduces_later_fuel_burn():
    from unicodes.aero import sizing

    mission = [("climb", 0.98), ("cruise", 0.9)]
    carried = sizing.estimate_weights(100.0, mission, payload=20)
    dropped = sizing.estimate_weights(100.0, mission, payload=20, payload_drop_after="climb")
    assert dropped.fuel_used == pytest.approx(2 + (98 - 20) * 0.1)
    assert dropped.fuel_used < carried.fuel_used


def test_take_off_weight_iteration_converges():
    from unicodes.aero import sizing

    mission = [("all", 0.7)]
    est = sizing.size_take_off_weight(mission, payload=10000, crew=700, A=0.2678, B=0.9979, W_guess=1e5)
    assert est.empty_tentative == pytest.approx(sizing.roskam_empty_weight(est.W_TO, 0.2678, 0.9979))


def test_airfoil_file_round_trip_and_pah(tmp_path):
    from unicodes.aero.airfoil import read_airfoil_file, relaxed_airfoil

    af = make_airfoil()
    path = tmp_path / "x.airfoil"
    af.write_airfoil_file(path)
    back = read_airfoil_file(path)
    assert back.alpha_c_l_max == pytest.approx(af.alpha_c_l_max)
    pah = relaxed_airfoil(af, 2.0)
    assert pah.c_l_alpha == pytest.approx(af.c_l_alpha / 2)
    assert pah.alpha_c_l_max == af.alpha_c_l_max
    assert pah.c_l(pah.alpha_c_l_max) == pytest.approx(af.c_l_max, abs=1e-6)


def test_camber_surface_slice_and_gradient():
    from unicodes.aero.airfoil import relaxed_airfoil
    from unicodes.aero.pah import CamberSurface, linear_camber_schedule

    lo = make_airfoil()
    hi = relaxed_airfoil(lo, 1.0)
    surf = CamberSurface.from_airfoils(lo, hi, 0.0, 1.0)
    # identical airfoils -> no camber dependence
    da, dc = surf.gradients("c_l")
    np.testing.assert_allclose(dc, 0, atol=1e-12)
    a = np.deg2rad(4)
    assert surf.slice("c_l", a, linear_camber_schedule(0, 0.5, 0)) == pytest.approx(lo.c_l(a), rel=1e-3)


POLAR_TEXT = """
       XFOIL         Version 6.99

 Calculated polar for: NACA 2412

 1 1 Reynolds number fixed          Mach number fixed

 xtrf =   1.000 (top)        1.000 (bottom)
 Mach =   0.000     Re =     0.500 e 6     Ncrit =   9.000

  alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr
 ------ -------- --------- --------- -------- -------- --------
  0.000   0.2442   0.00702   0.00175  -0.0533   0.6693   0.9969
  2.000   0.4697   0.00742   0.00191  -0.0539   0.5779   1.0000
"""


def test_xfoil_polar_parser_and_mesh():
    from unicodes.aero import xfoil

    p = xfoil.parse_polar(POLAR_TEXT)
    assert p.name == "NACA 2412" and p.Re == 5e5 and p.Ncrit == 9
    np.testing.assert_allclose(p.CL, [0.2442, 0.4697])

    def fake_run(name, alpha, Re):  # thin-airfoil stand-in for XFOIL
        m = int(name[4]) / 100
        return xfoil.Polar(name, Re, 9, alpha, 2 * np.pi * np.deg2rad(alpha) + 10 * m, *(np.zeros_like(alpha),) * 5)

    A, C, cl = xfoil.lift_mesh([0, 2], alpha_max=10, n_alpha=11, run=fake_run)
    assert A.shape == (21, 3)
    np.testing.assert_allclose(cl[:, 1], 2 * np.pi * A[:, 1], atol=1e-12)  # symmetric section in the middle
    np.testing.assert_allclose(cl[:, 0], -cl[::-1, 2])
