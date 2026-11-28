import numpy as np
import pytest

from unicodes.flight_dynamics.kinematics import lla_to_flat
from unicodes.fluids import friction_factor_colebrook, friction_factor_haaland, friction_factor_laminar
from unicodes.thermo import adiabatic_flame_temperature, build_relations, expand_stations, parse_formula, solve_relations
from unicodes.thermo.combustion import enthalpy_molar


def test_formula_parser():
    assert parse_formula("CH4O") == {"C": 1.0, "H": 4.0, "O": 1.0}
    assert parse_formula("C8.26H15.5") == {"C": 8.26, "H": 15.5}


def test_methane_thermo_and_flame_temperature():
    assert enthalpy_molar("CH4", 298.15) == pytest.approx(-74.6e3, rel=5e-3)
    assert enthalpy_molar("C3H8", 298.15) == pytest.approx(-103.9e3, rel=1e-2)
    T = adiabatic_flame_temperature("methane", 1.0)
    assert 2300 < T < 2350  # complete combustion, no dissociation
    assert adiabatic_flame_temperature("methane", 0.6) < T < adiabatic_flame_temperature("methane", 1.0, "O2")


def test_relation_network():
    eqs = ["F = m*a", "a = v/t"] + expand_stations(["pt_st/p_st = (1+(gamma-1)/2*M_st^2)^(gamma/(gamma-1))"], [0, 9])
    assert eqs[2].startswith("pt_0") and eqs[3].startswith("pt_9")
    out = solve_relations(build_relations(eqs), {"m": 2, "v": 10, "t": 5, "p_0": 1.0, "pt_0": 1.893, "gamma": 1.4})
    assert out["F"] == pytest.approx(4)
    assert out["M_0"] == pytest.approx(1.0, abs=1e-3)
    assert "M_9" not in out


def test_friction_factors():
    assert friction_factor_laminar(1000) == pytest.approx(0.064)
    for Re, e in ((1e4, 1e-5), (1e6, 1e-3), (1e7, 1e-2)):
        assert friction_factor_haaland(Re, e) == pytest.approx(friction_factor_colebrook(Re, e), rel=0.03)
    assert friction_factor_colebrook(1e5, 0) == pytest.approx(0.01799, rel=1e-3)  # smooth pipe


def test_lla_to_flat():
    n, e, d = lla_to_flat(38.96, -95.25, 300.0, 38.95, -95.26, -300.0)
    assert n == pytest.approx(1110.0, rel=5e-3)  # 0.01 deg of latitude
    assert e == pytest.approx(1110.0 * np.cos(np.deg2rad(38.95)), rel=5e-3)
    assert d == pytest.approx(0.0)
