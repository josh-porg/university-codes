import numpy as np
import pytest

from unicodes import propulsion_cycles as pc
from unicodes.propulsion_cycles import AIR, Gas


def ideal_turbojet_specific_thrust(M0, T0, pi_c, Tt4, g=1.4, cp=1004.0):
    """Mattingly's closed-form ideal turbojet (f << 1 neglected)."""
    R = cp * (g - 1) / g
    a0 = np.sqrt(g * R * T0)
    tau_r, tau_l, tau_c = 1 + (g - 1) / 2 * M0**2, Tt4 / T0, pi_c ** ((g - 1) / g)
    tau_t = 1 - tau_r / tau_l * (tau_c - 1)
    V9_a0 = np.sqrt(2 / (g - 1) * tau_l / (tau_r * tau_c) * (tau_r * tau_c * tau_t - 1))
    return a0 * (V9_a0 - M0)


def test_ideal_turbojet_matches_closed_form():
    # neglect the fuel mass: a huge heating value makes f -> 0
    r = pc.turbojet(2.0, 216.7, 19e3, 10, 1800, Qr=1e12, gas_c=AIR, gas_t=AIR)
    assert r.specific_thrust == pytest.approx(ideal_turbojet_specific_thrust(2.0, 216.7, 10, 1800), rel=1e-5)
    assert r.exits["9"].p == pytest.approx(19e3)


def test_turbojet_energy_and_efficiencies():
    r = pc.turbojet(0.8, 230, 30e3, 20, 1600, 42.8e6, pi_d=0.98, pi_b=0.95, pi_n=0.98, e_c=0.9, e_t=0.9,
                    eta_b=0.99, eta_m=0.99)
    # burner energy balance
    assert (1 + r.f) * 1156 * 1600 - 1004 * r.stations["3"].Tt == pytest.approx(0.99 * r.f * 42.8e6)
    assert 0 < r.eta_overall < r.eta_thermal < 1
    assert r.eta_overall == pytest.approx(r.specific_thrust * r.V0 / (r.f * 42.8e6))
    ab = pc.turbojet(0.8, 230, 30e3, 20, 1600, 42.8e6, Tt7=2000, gas_ab=Gas(1.3, 1240))
    assert ab.specific_thrust > r.specific_thrust and ab.tsfc > r.tsfc


def test_turbofan_reduces_to_turbojet():
    jet = pc.turbojet(0.8, 230, 30e3, 20, 1600, 42.8e6, e_c=0.9, e_t=0.9)
    fan = pc.turbofan(0.8, 230, 30e3, 20, 1.5, 0.0, 1600, 42.8e6, e_c=0.9, e_t=0.9)
    assert fan.specific_thrust == pytest.approx(jet.specific_thrust)


def test_mixer_of_identical_streams_is_lossless():
    s = pc.Station(200e3, 800.0)
    out, gas, M6A, M16 = pc.constant_area_mixer(1.0, s, AIR, 0.5, 2.0, s, AIR)
    assert out.pt == pytest.approx(200e3) and out.Tt == pytest.approx(800)
    assert M6A == pytest.approx(0.5) and M16 == pytest.approx(0.5)


def test_mixed_turbofan_matches_pressures_and_conserves_energy():
    r = pc.mixed_turbofan(1.6, 216.7, 10e3, 20, 3.0, 1700, 42.8e6, 0.5, e_c=0.9, e_f=0.89, e_t=0.89,
                          pi_b=0.96, eta_b=0.99, eta_m=0.99)
    st = r.stations
    assert st["16"].pt == pytest.approx(st["5"].pt)
    assert r.alpha > 0
    a, f = r.alpha, r.f
    lhs = (1 + f) * 1156 * st["5"].Tt + a * 1004 * st["13"].Tt
    assert lhs == pytest.approx((1 + f + a) * r.extra["cp_6A"] * st["6A"].Tt)
    assert r.extra["pi_M_ideal"] < 1.0


def test_convergent_nozzle_chokes():
    e = pc.nozzle(3.0e5, 800, 1e5, AIR, convergent=True)
    assert e.choked and e.M == pytest.approx(1.0) and e.p > 1e5
    e = pc.nozzle(1.5e5, 800, 1e5, AIR, convergent=True)
    assert not e.choked and e.p == 1e5 and e.M < 1


def test_efficiency_inversions():
    tau, eta = pc.compressor(20, 0.9)
    eta2, e = pc.compressor_efficiencies(20, tau)
    assert eta2 == pytest.approx(eta) and e == pytest.approx(0.9)
    pi, eta_t = pc.turbine(0.7, 0.88)
    assert pc.turbine_efficiencies(pi, 0.7) == pytest.approx((eta_t, 0.88))


def test_compressor_stage_ideal():
    s = pc.compressor_stage(0.4, 774.926, 150, 0.75, 287, 100e3, 1.5)
    assert s.ct2 == pytest.approx(2 * s.U * 0.25)
    assert s.efficiency == pytest.approx(1.0)
    lossy = pc.compressor_stage(0.4, 774.926, 150, 0.75, 287, 100e3, 1.5, 0.05, 0.05)
    assert lossy.efficiency < 1 and lossy.Tt3 == pytest.approx(s.Tt3)
