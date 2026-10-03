import numpy as np
import pytest

from unicodes.composites import Laminate, Ply, mechanical_properties

GRAPHITE_EPOXY = dict(E_1=181e9, E_2=10.3e9, G_12=7.17e9, nu_12=0.28)  # T300/5208, Kaw


def test_reduced_stiffness_matches_textbook():
    # Kaw, Mechanics of Composite Materials, Example 2.6
    Q = Ply(**GRAPHITE_EPOXY).Q
    np.testing.assert_allclose(
        Q / 1e9, [[181.8, 2.897, 0], [2.897, 10.35, 0], [0, 0, 7.17]], rtol=1e-3, atol=1e-9
    )


def test_qbar_at_zero_is_q():
    ply = Ply(**GRAPHITE_EPOXY)
    np.testing.assert_allclose(ply.Qbar(0.0), ply.Q, atol=1e-3)


def test_isotropic_laminate_recovers_material_properties():
    E, nu = 70e9, 0.3
    ply = Ply(E, E, E / (2 * (1 + nu)), nu)
    lam = Laminate([ply], None, 1e-3, np.deg2rad([0, 45, -45, 90]), symmetry=2)
    assert lam.E_x == pytest.approx(E)
    assert lam.nu_xy == pytest.approx(nu)
    assert lam.t_laminate == pytest.approx(8e-3)


def test_symmetric_laminate_has_no_coupling():
    lam = Laminate([Ply(**GRAPHITE_EPOXY)], None, 0.125e-3, np.deg2rad([0, 30, -45]), symmetry=2)
    np.testing.assert_allclose(lam.B, 0, atol=1e-6 * np.abs(lam.A).max())
    np.testing.assert_allclose(lam.ply_angles, np.deg2rad([0, 30, -45, -45, 30, 0]))


def test_antisymmetric_mirror_flips_angles():
    lam = Laminate([Ply(**GRAPHITE_EPOXY)], None, 1e-4, np.deg2rad([45, 90]), symmetry=-2)
    np.testing.assert_allclose(lam.ply_angles, np.deg2rad([45, 90, 90, -45]))


def test_degrees_are_rejected():
    with pytest.raises(ValueError):
        Laminate([Ply(**GRAPHITE_EPOXY)], None, 1e-4, [0, 90], symmetry=0)


def test_first_ply_failure_under_axial_tension():
    ply = Ply(**GRAPHITE_EPOXY, strain_allowables=[0.01, 0.005, -0.01, -0.02, 0.02])
    lam = Laminate([ply], None, 0.125e-3, np.deg2rad([0, 90]), symmetry=2)
    loads = np.array([4e5, 0, 0, 0, 0, 0])  # N/m, ~0.8% axial strain
    failed, failures, ratios = lam.load(loads, "Max_Strain")
    # 90-degree plies see the axial strain as transverse strain, which has the lower allowable
    assert failed.tolist() == [False, True, True, False]
    assert {f.failure_mode for f in failures} == {2}
    assert "transverse tensile strain" in str(failures[0])


def test_mechanical_properties_rule_of_mixtures():
    E_1, E_2, nu_12, G_12 = mechanical_properties(230e9, 22e9, 0.2, 3.5e9, 1.3e9, 0.35, 0.6, 2, 1)
    assert E_1 == pytest.approx(0.6 * 230e9 + 0.4 * 3.5e9)
    assert nu_12 == pytest.approx(0.26)
    assert 3.5e9 < E_2 < 230e9
