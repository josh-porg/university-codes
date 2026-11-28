import numpy as np
import pytest
from scipy.integrate import trapezoid

from unicodes import atmosphere, gasdynamics as gd
from unicodes.aero import boundary_layer as bl
from unicodes.aero import wing
from unicodes.io import to_latex_table


def test_isa_reference_values():
    s = atmosphere.isa([0, 11000, 20000, 32000])
    np.testing.assert_allclose(s.T, [288.15, 216.65, 216.65, 228.65], atol=1e-9)
    np.testing.assert_allclose(s.p, [101325, 22632.06, 5474.89, 868.02], rtol=2e-5)
    assert s.rho[0] == pytest.approx(1.225, rel=1e-4)
    assert s.a[0] == pytest.approx(340.294, rel=1e-4)
    assert s.mu[0] == pytest.approx(1.7894e-5, rel=1e-3)


def test_isa_hot_day_keeps_pressure():
    std, hot = atmosphere.isa(3000), atmosphere.isa(3000, delta_T=10)
    assert hot.p == pytest.approx(std.p)
    assert hot.rho == pytest.approx(std.rho * std.T / hot.T)


def test_isentropic_and_area_ratio():
    r = gd.isentropic(2.0)
    assert r.p0_p == pytest.approx(7.8244, rel=1e-4)
    assert r.T0_T == pytest.approx(1.8)
    assert r.A_Astar == pytest.approx(1.6875, rel=1e-4)
    assert gd.mach_from_area_ratio(1.6875, supersonic=True) == pytest.approx(2.0, rel=1e-6)
    assert gd.mach_from_area_ratio(1.6875, supersonic=False) == pytest.approx(0.3722, rel=1e-3)


def test_normal_shock_m2():
    s = gd.normal_shock(2.0)
    assert s.M2 == pytest.approx(0.5774, rel=1e-4)
    assert s.p2_p1 == pytest.approx(4.5)
    assert s.rho2_rho1 == pytest.approx(2.6667, rel=1e-4)
    assert s.p02_p01 == pytest.approx(0.7209, rel=1e-4)
    assert gd.mach_from_pitot_supersonic(gd.pitot_rayleigh(2.0), 1.0) == pytest.approx(2.0)


def test_oblique_shock_textbook_case():
    beta = gd.wave_angle(3.0, np.deg2rad(20))
    assert np.rad2deg(beta) == pytest.approx(37.76, abs=0.05)
    s = gd.oblique_shock(3.0, beta, np.deg2rad(20))
    assert s.M2 == pytest.approx(1.99, abs=0.01)
    assert np.isnan(gd.wave_angle(1.5, np.deg2rad(30)))


def test_prandtl_meyer_round_trip():
    assert np.rad2deg(gd.prandtl_meyer(2.0)) == pytest.approx(26.38, abs=0.01)
    assert gd.expansion_fan(2.0, np.deg2rad(10)) == pytest.approx(2.385, abs=2e-3)


def test_critical_cp_and_mcrit():
    assert gd.critical_pressure_coefficient(1.0) == pytest.approx(0, abs=1e-12)
    M = wing.critical_mach_from_cp(-0.43, "prandtl_glauert")
    assert wing.prandtl_glauert(-0.43, M) == pytest.approx(gd.critical_pressure_coefficient(M))


def test_lift_slopes_agree_in_limits():
    A = 8.0
    assert wing.helmbold_lift_slope(A) == pytest.approx(wing.rectangular_wing_lift_slope(A))
    a0 = 2 * np.pi
    assert wing.finite_wing_lift_slope(a0, 1e9) == pytest.approx(a0)
    assert wing.rescale_lift_slope(4.0, 6.0, 6.0) == pytest.approx(4.0)
    assert wing.polhamus_lift_slope(a0, A, 0.0, 1.0, 0.5) > wing.polhamus_lift_slope(a0, A, 0.0, 1.0)


def test_mgc_and_schrenk_integrates_to_half_lift():
    c_bar, _, y = wing.mean_geometric_chord(2.0, 10.0, 1.0)
    assert c_bar == pytest.approx(2.0) and y == pytest.approx(2.5)
    assert wing.root_chord_from_mgc(c_bar, 1.0) == pytest.approx(2.0)
    y = np.linspace(0, 5, 20001)
    L = wing.schrenk_lift_distribution(y, 10, 2.0, 0.8, 1000.0)
    assert trapezoid(L, y) == pytest.approx(500, rel=1e-3)


def test_flat_plate_drag_laminar_limit():
    nu, rho, V = 1.5e-5, 1.2, 10.0
    D = bl.flat_plate_drag(V, 1.0, 1.0, rho, nu, Re_crit=1e9)
    assert D == pytest.approx(2 * 0.5 * rho * V**2 * bl.cf_laminar(V / nu))
    assert bl.flat_plate_drag(V, 1.0, 1.0, rho, nu, Re_crit=5e5) > D


def test_latex_table():
    text = to_latex_table([[1, 2.5], [float("inf"), 3]], headers=["a", "b"])
    assert r"\infty" in text and text.startswith(r"\begin{tabular}{ll}")


from unicodes.aero import naca, pressure, thin_airfoil  # noqa: E402


def test_naca_shapes():
    af = naca.naca4("2412", n=60)
    x_u, z_u = af.upper
    assert x_u[0] == pytest.approx(1) and x_u[-1] == pytest.approx(0)
    t = 0.12
    assert (af.z.max() - af.z.min()) == pytest.approx(t, rel=0.1)
    sym = naca.naca4("0012", n=40)
    np.testing.assert_allclose(sym.upper[1][::-1], -sym.lower[1], atol=1e-12)
    a5 = naca.naca5("23012", n=24)
    assert len(a5.x) == 49
    assert a5.z_camber.max() == pytest.approx(0.0184, abs=1e-3)


def test_flat_plate_pressure_integration():
    x = np.linspace(0, 1, 101)
    upper, lower = (x, 0 * x), (x, 0 * x)
    res = pressure.section_coefficients(0.0, -np.ones(100), np.ones(100), upper, lower)
    assert res.c_l == pytest.approx(2.0)
    assert res.c_m_le == pytest.approx(-1.0)
    assert res.c_m_qc == pytest.approx(-0.5)


def test_thin_airfoil_symmetric_and_cambered():
    flat = thin_airfoil.ThinAirfoil(lambda x: 0 * x)
    assert flat.zero_lift_angle == pytest.approx(0, abs=1e-12)
    assert flat.c_l(0.1) == pytest.approx(2 * np.pi * 0.1)
    # NACA 2412 camber: alpha_L0 about -2.08 deg, cm_c/4 about -0.053
    m, p = 0.02, 0.4
    slope = lambda x: np.where(x < p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))  # noqa: E731
    af = thin_airfoil.ThinAirfoil(slope, (p,))
    assert np.rad2deg(af.zero_lift_angle) == pytest.approx(-2.08, abs=0.01)
    assert af.c_m_quarter_chord == pytest.approx(-0.053, abs=1e-3)
