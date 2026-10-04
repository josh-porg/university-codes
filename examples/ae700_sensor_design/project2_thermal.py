"""AE 700 project 2: thermal-infrared line scanner sizing (``Project_2_algegbra``).

Find the detection wavelength at which a diffraction-limited detector
element of a 3000-element array produces the minimum signal, then size the
aperture for the required ground resolution. The symbolic ``vpasolve``
becomes ``brentq``.
"""

import numpy as np
from scipy.optimize import brentq

from unicodes import remote_sensing as rs
from unicodes.units import fahrenheit_to_kelvin

T = fahrenheit_to_kelvin(53.33)
epsilon, n_det, responsivity, tau_opt, S_min = 0.7, 3000, 1e9, 0.96, 5e5
SW_min = 6 * 1609.344
R_max = 3.5
R_min = SW_min / n_det
for label, R in (("min", R_min), ("mid", (R_min + R_max) / 2), ("max", R_max)):
    H, f = 4850.0, 0.5
    SW = R * n_det
    lam_peak = rs.wien_peak_wavelength(T)

    # S = Resp * L * (pi/4) d^2 cos^4(R/H) tau / f^2 * (1.22 f lambda / d)^2, d cancels:
    def signal(lam):
        return responsivity * rs.planck_spectral_radiance(T, lam, epsilon) * 1.22**2 * lam**2 * np.pi / 4 \
            * np.cos(R / H) ** 4 * tau_opt

    lam_det = brentq(lambda lam: signal(lam) - S_min, 3e-6, lam_peak)
    d = 1.22 * H * lam_det / R
    dx = 1.22 * f / d * lam_det
    print(f"{label} resolution {R:.3f} m (swath {SW:.0f} m): lambda_det = {lam_det * 1e6:.3f} um "
          f"({(lam_peak - lam_det) / lam_peak:.1%} below Wien peak {lam_peak * 1e6:.2f} um), aperture {d:.4f} m, "
          f"detector pitch {dx * 1e6:.2f} um, array width {n_det * dx * 1e3:.1f} mm")
