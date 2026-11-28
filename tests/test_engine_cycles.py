import numpy as np
import pytest

from unicodes.thermo import enthalpy_molar
from unicodes.thermo.engine_cycles import dual_cycle, dual_cycle_pv, dual_cycle_work, heat_temperature


def test_dual_cycle_work_matches_the_pv_loop_area():
    p1, v1, eps, alpha, beta, n1, n2 = 100e3, 0.85, 16.0, 1.2, 2.5, 1.38, 1.25
    v, p = dual_cycle_pv(p1, v1, eps, alpha, beta, n1, n2, n=20001)
    area = np.trapezoid(p, v) if hasattr(np, "trapezoid") else np.trapz(p, v)  # clockwise loop: integral of p dv
    assert dual_cycle_work(p1, v1, eps, alpha, beta, n1, n2) == pytest.approx(area, rel=1e-6)


def test_heat_temperature_closes_the_energy_balance():
    moles = {"N2": 3.76, "O2": 1.0}
    T = heat_temperature(moles, 800.0, 1e5, "v")
    du = sum(n * (enthalpy_molar(s, T) - enthalpy_molar(s, 800.0)) for s, n in moles.items()) - 4.76 * 8.314 * (T - 800)
    assert du == pytest.approx(1e5)
    assert heat_temperature(moles, 800.0, 1e5, "p") < T


def test_dual_cycle_project_engine():
    r = dual_cycle("gasoline (heavy)")
    assert r.lhv_fuel == pytest.approx(43.7e6, rel=0.01)
    assert r.alpha == pytest.approx(5.5e6 / (100e3 * 16**1.38))
    assert 0.4 < r.efficiency < 0.6 and r.beta > 1
    otto = dual_cycle("gasoline (heavy)", p_max=50e6)
    assert otto.beta == 1.0 and otto.p[3] < 50e6
