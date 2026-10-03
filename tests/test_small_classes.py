import math

import pytest

from unicodes.fluids import PipeSection
from unicodes.materials import Material
from unicodes.propulsion import RocketBullet


def test_rocket_bullet_derived_values():
    rb = RocketBullet("test", 0.0127, 0.006, 0.06, 0.02, 0.002, 0.011, 0.05, 200, 0.01, 0.0, 0.5)
    assert rb.nozzle_exit_area == pytest.approx(math.pi * 0.003**2)
    assert rb.dry_mass == pytest.approx(0.04)
    assert rb.burn_time == 0.5


def test_pipe_section_continuity_and_loss():
    pipe = PipeSection(area_initial=2.0, area_final=1.0, K_L=0.5, rho=1000, V_initial=1.0)
    assert pipe.V_final == 2.0
    assert pipe.pressure_loss == pytest.approx(250)
    assert pipe.delta_p == pytest.approx(1000 * (1 - 4) / 2 - 250)
    with pytest.raises(ValueError):
        PipeSection(1, 1, 0.1, 1000)


def test_material_specific_properties():
    al = Material("6061-T6", 68.9e9, 310e6, 276e6, 2700, 5.0)
    assert al.specific_stiffness == pytest.approx(68.9e9 / 2700)
