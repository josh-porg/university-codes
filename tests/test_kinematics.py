import numpy as np
import pytest

from unicodes.flight_dynamics import kinematics as k


def test_quaternion_matches_dcm():
    angles = np.deg2rad([10, 5, -45])
    H = k.dcm_inertial_to_body(angles)
    q = k.quat_from_euler(*angles)
    v = np.array([850.0, -30.0, 20.0])
    np.testing.assert_allclose(k.rotate_active(v, q), H @ v, atol=1e-9)
    np.testing.assert_allclose(k.rotate_passive(H @ v, q), v, atol=1e-9)
    np.testing.assert_allclose(k.euler_from_quat(q), angles, atol=1e-12)
    np.testing.assert_allclose(H @ H.T, np.eye(3), atol=1e-12)


def test_axis_angle_and_inverse():
    q = k.quat_from_axis_angle([0, 0, 1], np.pi / 2)
    np.testing.assert_allclose(k.quat_multiply(q, k.quat_inverse(q)), [1, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(k.quat_from_euler(0, 0, np.pi / 2), q, atol=1e-12)


def test_path_and_flow_angles():
    V = k.velocity_from_path_angles(450, np.deg2rad(50), np.deg2rad(25))
    g, chi = k.flight_path_angles(V)
    assert np.rad2deg(g) == pytest.approx(50) and np.rad2deg(chi) == pytest.approx(25)
    a, b = k.air_flow_angles([100, 0, 10])
    assert a == pytest.approx(np.arctan(0.1)) and b == 0


def test_aerodynamic_centre_wing_only():
    assert k.aircraft_aerodynamic_centre(0.25, 5, 4.0, 4, 0.0, 1, 0.3, 10) == pytest.approx(0.25)
