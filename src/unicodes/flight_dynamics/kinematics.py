"""Attitude kinematics: direction cosines, quaternions, flight-path and air-flow angles (AE 550).

Quaternions are ``[q0, q1, q2, q3]`` scalar first. Euler angles are the
aerospace 3-2-1 sequence ``(phi, theta, psi)`` = roll, pitch, yaw. The
MATLAB ``eulerAngles2quaternion``/``quaternion2eulerAngles`` pair had roll
and pitch swapped and ``axisAngle2qaternion`` swapped sin and cos; these
follow the standard definitions, and the multiplication is the Hamilton
product (see :func:`quat_multiply`).
"""

from __future__ import annotations

import numpy as np


def dcm_inertial_to_body(phi, theta=None, psi=None):
    """Direction cosine matrix H_IB taking inertial (NED) components to body components (``HIB``).

    Accepts ``(phi, theta, psi)`` or a single 3-vector.
    """
    if theta is None:
        phi, theta, psi = phi
    cps, sps = np.cos(psi), np.sin(psi)
    cth, sth = np.cos(theta), np.sin(theta)
    cph, sph = np.cos(phi), np.sin(phi)
    H_I1 = np.array([[cps, sps, 0], [-sps, cps, 0], [0, 0, 1]])
    H_12 = np.array([[cth, 0, -sth], [0, 1, 0], [sth, 0, cth]])
    H_2B = np.array([[1, 0, 0], [0, cph, sph], [0, -sph, cph]])
    return H_2B @ H_12 @ H_I1


def quat_multiply(r, s):
    """Hamilton product ``r * s`` (``quatMult``).

    The MATLAB ``quatMult`` flipped the sign of the cross-product term, so it
    computed ``s * r``.
    """
    r0, r1, r2, r3 = r
    s0, s1, s2, s3 = s
    return np.array([
        r0 * s0 - r1 * s1 - r2 * s2 - r3 * s3,
        r0 * s1 + r1 * s0 + r2 * s3 - r3 * s2,
        r0 * s2 - r1 * s3 + r2 * s0 + r3 * s1,
        r0 * s3 + r1 * s2 - r2 * s1 + r3 * s0,
    ])  # fmt: skip


def quat_conjugate(q):
    return np.array([q[0], -q[1], -q[2], -q[3]], dtype=float)


def quat_inverse(q):
    return quat_conjugate(q) / np.dot(q, q)


def quat_from_axis_angle(axis, angle):
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    return np.concatenate([[np.cos(angle / 2)], np.sin(angle / 2) * axis])


def quat_from_euler(phi, theta, psi):
    """Quaternion of the 3-2-1 rotation inertial -> body."""
    cph, sph = np.cos(phi / 2), np.sin(phi / 2)
    cth, sth = np.cos(theta / 2), np.sin(theta / 2)
    cps, sps = np.cos(psi / 2), np.sin(psi / 2)
    return np.array([
        cph * cth * cps + sph * sth * sps,
        sph * cth * cps - cph * sth * sps,
        cph * sth * cps + sph * cth * sps,
        cph * cth * sps - sph * sth * cps,
    ])  # fmt: skip


def euler_from_quat(q):
    """``(phi, theta, psi)`` from a unit quaternion; handles the pitch = +-90 deg singularity."""
    q0, q1, q2, q3 = q
    s = np.clip(2 * (q0 * q2 - q1 * q3), -1, 1)
    theta = np.arcsin(s)
    if np.isclose(abs(s), 1):
        return 0.0, theta, -2 * np.sign(s) * np.arctan2(q1, q0)
    phi = np.arctan2(2 * (q0 * q1 + q2 * q3), q0**2 - q1**2 - q2**2 + q3**2)
    psi = np.arctan2(2 * (q0 * q3 + q1 * q2), q0**2 + q1**2 - q2**2 - q3**2)
    return phi, theta, psi


def rotate_active(v, q):
    """Rotate vector ``v`` by quaternion ``q``: ``q* (0, v) q`` (``rotateVectorWithQuaternion``/``quatRot``).

    With ``q = quat_from_euler(...)`` this gives the body components of an
    inertial vector, i.e. ``dcm_inertial_to_body(...) @ v``.
    """
    p = np.concatenate([[0.0], v])
    return quat_multiply(quat_multiply(quat_conjugate(q), p), q)[1:]


def rotate_passive(v, q):
    """``q (0, v) q*`` (``quatRotPassive``), the inverse of :func:`rotate_active`."""
    p = np.concatenate([[0.0], v])
    return quat_multiply(quat_multiply(q, p), quat_conjugate(q))[1:]


def air_flow_angles(V_body):
    """Angle of attack and sideslip ``(alpha, beta)`` from body-axis air velocity (``getAirFlowAngles``)."""
    u, v, w = V_body
    return np.arctan2(w, u), np.arcsin(v / np.linalg.norm(V_body))


def flight_path_angles(V_inertial):
    """Climb angle and heading ``(gamma, chi)`` from NED velocity.

    (``getFlightPathAngles`` returned ``atan2(v_y, v_x)`` as gamma and the
    climb angle as xi, i.e. the names swapped.)
    """
    vx, vy, vz = V_inertial
    return np.arctan2(-vz, np.hypot(vx, vy)), np.arctan2(vy, vx)


def velocity_from_path_angles(V, gamma, chi):
    """NED velocity from speed, climb angle and heading (AE 550 HW 2 problem 5)."""
    return V * np.array([np.cos(gamma) * np.cos(chi), np.cos(gamma) * np.sin(chi), -np.sin(gamma)])


def aircraft_aerodynamic_centre(x_ac_wf, C_L_alpha_wf, x_ac_h, C_L_alpha_h, S_h, eta_h, deda_h, S,
                                x_ac_c=0.0, C_L_alpha_c=0.0, S_c=0.0, eta_c=1.0, deda_c=0.0, munk_shift=0.0):
    """Whole-aircraft a.c. from wing-fuselage, tail and canard contributions (``aircraft_aerodynamic_center_derivation``)."""
    wf = C_L_alpha_wf
    h = eta_h * S_h / S * C_L_alpha_h * (1 - deda_h)
    c = eta_c * S_c / S * C_L_alpha_c * (1 + deda_c)
    return (wf * (x_ac_wf + munk_shift) + h * x_ac_h + c * x_ac_c) / (wf + h + c)
