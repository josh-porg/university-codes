"""Unit conversion factors (replaces ``convlength``, ``convforce``, ``symunit`` ...).

Multiply a value in the named unit by the constant to get SI::

    from unicodes.units import FT, LBF
    span_m = 110 * FT
    weight_lbf = weight_N / LBF
"""

FT = 0.3048  # m
INCH = 0.0254  # m
MILE = 1609.344  # m
NMI = 1852.0  # m
LBF = 4.4482216152605  # N
LBM = 0.45359237  # kg
SLUG = LBF / FT  # kg (14.5939...)
SLUG_PER_FT3 = SLUG / FT**3  # kg/m^3
PSF = LBF / FT**2  # Pa
PSI = LBF / INCH**2  # Pa
KT = NMI / 3600  # m/s
MPH = MILE / 3600  # m/s
KMH = 1000 / 3600  # m/s
FPM = FT / 60  # m/s (ft/min)
HP = 550 * FT * LBF  # W
RANKINE = 5 / 9  # K per degree R
DEG = 3.141592653589793 / 180  # rad
LBF_S_PER_FT2 = LBF / FT**2  # Pa s (viscosity)
G0 = 9.80665  # m/s^2


def fahrenheit_to_kelvin(T_F):
    return (T_F - 32) * 5 / 9 + 273.15


def celsius_to_kelvin(T_C):
    return T_C + 273.15
