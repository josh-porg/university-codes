import pytest

from unicodes.thermo import GasLaw, Relation, atom_balance, cp_molar, enthalpy_molar, heat_of_combustion


@pytest.mark.parametrize(
    "species, h_formation",  # J/mol at 298 K
    [("CO2", -393.5e3), ("H2O", -241.8e3), ("O2", 0.0), ("N2", 0.0), ("H2", 0.0)],
)
def test_nasa_enthalpy_recovers_heat_of_formation(species, h_formation):
    assert enthalpy_molar(species, 298.15) == pytest.approx(h_formation, abs=0.5e3)


def test_nasa_polynomials_continuous_at_1000K():
    for s in ("CO2", "H2O", "O2", "N2"):
        assert enthalpy_molar(s, 999.999) == pytest.approx(enthalpy_molar(s, 1000.0), rel=1e-3)
        assert cp_molar(s, 999.999) == pytest.approx(cp_molar(s, 1000.0), rel=1e-2)


def test_atom_balance_lean_gasoline():
    rxn = atom_balance("Gasoline (Heavy)", 0.8)
    assert rxn.a == pytest.approx((8.26 + 15.5 / 4) / 0.8)
    # oxygen atoms balance
    assert 2 * rxn.a == pytest.approx(2 * rxn.b + rxn.c + 2 * rxn.d)


def test_gasoline_lower_heating_value():
    rxn = atom_balance("Gasoline (Heavy)", 1.0)
    lhv_mj_per_kg = heat_of_combustion(rxn, 298.0) / rxn.fuel.molar_mass / 1e3  # J/g = kJ/kg
    assert 42 < lhv_mj_per_kg < 46  # ~44 MJ/kg


def test_hydrogen_lower_heating_value():
    rxn = atom_balance("hydrogen", 1.0)
    assert heat_of_combustion(rxn, 298.0) == pytest.approx(241.8e3, rel=5e-3)


def test_ideal_gas_relations():
    import sympy as sp

    gl = GasLaw("p*v = R*T", cp="cp0")
    R, cp0 = sp.symbols("R cp0")
    assert sp.simplify(gl.partials["PrTcV"] - R / gl.v) == 0
    assert gl.in_terms_of(gl.cv) == cp0 - R
    assert gl.in_terms_of(gl.partials["TrPcH"]) == 0  # no Joule-Thomson effect
    assert gl.reciprocity_ok


def test_van_der_waals_internal_pressure():
    import sympy as sp

    gl = GasLaw("(p + a/v**2)*(v - b) = R*T", cp="cp0")
    a = sp.Symbol("a")
    assert sp.simplify(gl.in_terms_of(gl.partials["UrVcT"]) - a / gl.v**2) == 0


def test_relation_solves_for_missing_variable():
    r = Relation("p*v = R*T", positive=["p", "v", "R", "T"])
    (v,) = r.solve(p=101325, R=287, T=288)
    assert float(v) == pytest.approx(287 * 288 / 101325)
