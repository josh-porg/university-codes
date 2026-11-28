"""Optical, thermal-infrared and SAR sensor sizing relations (AE 700 vehicle design requirements).

Ports the helper functions of ``project_1_v_1_0``/``Project1Runthrough_5``
(``res``, ``Swath``, ``ReflecRad``, ``falloff``, ``signal``, ``ToD2theta``),
``lambda_max_blackbody``, ``spectral_radiance_labertian_blackbody`` and the
radar relations of ``Project_3_algrebra``/``final_exam_problem_2``. SI units.
"""

from __future__ import annotations

import numpy as np

H_PLANCK = 6.62607015e-34  # the MATLAB used 6.625e-34
K_BOLTZMANN = 1.380649e-23
C_LIGHT = 299792458.0
WIEN = 2.898e-3  # m K


def wien_peak_wavelength(T):
    """Wavelength of peak blackbody emission ``b / T`` (``lambda_max_blackbody``)."""
    return WIEN / np.asarray(T, dtype=float)


def planck_spectral_radiance(T, wavelength, emissivity=1.0):
    """Spectral radiance of a Lambertian grey body, W/(m^2 sr m) (``spectral_radiance_labertian_blackbody``)."""
    lam = np.asarray(wavelength, dtype=float)
    return emissivity * 2 * H_PLANCK * C_LIGHT**2 / lam**5 / np.expm1(H_PLANCK * C_LIGHT / (lam * K_BOLTZMANN * T))


def rayleigh_angle(wavelength, aperture):
    """Diffraction-limited angular resolution ``1.22 lambda / D`` (rad)."""
    return 1.22 * wavelength / aperture


def ground_resolution(altitude, wavelength, aperture):
    """Diffraction-limited ground resolution ``1.22 H lambda / D``."""
    return altitude * rayleigh_angle(wavelength, aperture)


def altitude_for_resolution(resolution, wavelength, aperture):
    """Altitude giving ``resolution`` with aperture ``D`` (``res`` returned H/D)."""
    return resolution * aperture / (1.22 * wavelength)


def swath_width(altitude, detector_width, focal_length):
    """Push-broom swath width ``H w_d / f`` (``Swath``)."""
    return altitude * detector_width / focal_length


def sun_zenith_from_time(hour):
    """Sun angle from local noon, degrees, for a decimal 24 h time (``ToD2theta``)."""
    return hour * 15.0 - 180.0


def reflected_radiance(solar_irradiance, reflectance, sun_angle_deg):
    """Radiance reflected by a Lambertian surface ``E r cos(theta) / pi`` (``ReflecRad``)."""
    return solar_irradiance * reflectance * np.cos(np.deg2rad(sun_angle_deg)) / np.pi


def image_irradiance(radiance, aperture_area, focal_length, off_axis_angle, transmission):
    """Focal-plane irradiance with cos^4 fall-off ``L A cos^4(theta) tau / f^2`` (``falloff``).

    ``off_axis_angle`` in radians; the MATLAB passed the Rayleigh angle (rad) to ``cosd``.
    """
    return radiance * aperture_area * np.cos(off_axis_angle) ** 4 / focal_length**2 * transmission


def detector_signal(responsivity, irradiance, detector_area):
    """Detector output ``R E A_det`` (``signal``)."""
    return responsivity * irradiance * detector_area


# Synthetic aperture radar


def range_resolution(pulse_width, look_angle):
    """Ground-range resolution ``c tau / (2 sin theta)``."""
    return C_LIGHT * pulse_width / (2 * np.sin(look_angle))


def azimuth_resolution(antenna_length):
    """Focused SAR azimuth resolution ``D / 2``."""
    return antenna_length / 2


def minimum_prf(velocity, antenna_length):
    """Minimum pulse repetition frequency ``2 V / D``."""
    return 2 * velocity / antenna_length


def radar_received_power(P_avg, G_t, G_r, wavelength, rcs, slant_range, pulse_width, prf):
    """Radar equation with average power: ``P_avg G_t G_r lambda^2 sigma / ((4 pi)^3 rho^4 tau PRF)``."""
    return P_avg * G_t * G_r * wavelength**2 * rcs / ((4 * np.pi) ** 3 * slant_range**4 * pulse_width * prf)


def rcs_for_received_power(P_r, P_avg, G_t, G_r, wavelength, slant_range, pulse_width, prf):
    """Radar cross-section giving received power ``P_r`` (inverse of :func:`radar_received_power`)."""
    return P_r * (4 * np.pi) ** 3 * slant_range**4 * pulse_width * prf / (P_avg * G_t * G_r * wavelength**2)


def resolution_cell_rcs(range_res, azimuth_res, look_angle, wavelength):
    """RCS of a flat resolution cell ``4 pi R_r^2 R_a^2 cos(theta) / lambda^2`` (Project 3)."""
    return 4 * np.pi * range_res**2 * azimuth_res**2 * np.cos(look_angle) / wavelength**2
