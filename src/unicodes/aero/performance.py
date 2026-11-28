"""Roskam Part I performance constraints, climb and power (AE 521, AE 722).

All inputs and outputs are SI (N, m, s, W, N/m^2); where Roskam's
correlations are in British units the conversion happens inside. Functions
work elementwise so constraint diagrams can be drawn from arrays of wing
loading.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ..atmosphere import RHO0, isa
from ..units import FPM, FT, HP, LBF, PSF, SLUG_PER_FT3
from .wing import max_lift_to_drag, parabolic_drag_polar

# Mission fuel fractions and weights (Roskam Part I ch. 2)


def breguet_range(V, L_D, I_sp, W_initial, W_final):
    """Breguet range ``V (L/D) I_sp ln(W_i/W_f)`` with I_sp in seconds of weight flow (``BreguetRange``)."""
    return V * L_D * I_sp * np.log(W_initial / W_final)


def breguet_endurance(L_D, I_sp, W_initial, W_final):
    """Breguet endurance ``(L/D) I_sp ln(W_i/W_f)`` (``BreguetEndurance``)."""
    return L_D * I_sp * np.log(W_initial / W_final)


def fuel_fraction_range(distance, V, I_sp, L_D):
    """W_end/W_start of a cruise leg from the Breguet range equation (``FuelFraction_Range``)."""
    return np.exp(-distance / (V * L_D * I_sp))


def fuel_fraction_endurance(duration, I_sp, L_D):
    """W_end/W_start of a loiter leg (``FuelFraction_Endurance``)."""
    return np.exp(-duration / (L_D * I_sp))


def crew_weight(n_crew, military: bool = False):
    """Crew weight (N): 200 lbf military, 175 + 30 lbf baggage civil (``Weight_Crew``)."""
    return n_crew * (200 if military else 205) * LBF


def passenger_weight(n_passengers, long_distance: bool = False):
    """Passengers plus baggage (N): 175 lbf + 40/30 lbf bags (``Weight_PassengerAndBaggage``)."""
    return n_passengers * (175 + (40 if long_distance else 30)) * LBF


def empty_weight(W_TO, A, B):
    """Roskam statistical empty weight (N) from regression constants A, B (``Weight_Empty``).

    The course read A and B from ``Empty Weight Correlation Coefficients.csv``
    (Roskam Part I Table 2.15); pass the row for the airplane type.
    """
    return 10 ** ((np.log10(np.asarray(W_TO) / LBF) - A) / B) * LBF


# Roskam Part I Table 3.5 wetted-area regression ``log10 S_wet = c + d log10 W_TO`` (ft^2, lbf)
WETTED_AREA_COEFFS = {
    "Homebuilts": (1.2362, 0.4319),
    "Single Engine Propeller Driven": (1.0892, 0.5147),
    "Twin Engine Propeller Driven": (0.8635, 0.5632),
    "Agricultural": (1.0477, 0.5326),
    "Business Jets": (0.2263, 0.6977),
    "Regional Turbojets": (-0.0866, 0.8099),
    "Transport Jets": (0.0199, 0.7531),
    "Military Trainers": (0.8565, 0.5423),
    "Fighters": (-0.1289, 0.7506),
    "Mil. Patrol, Bomb and Transport": (0.1628, 0.7316),
    "Flying Boats, Amph. and Float": (0.6295, 0.6708),
    "Supersonic Cruise Airplanes": (-1.1868, 0.9609),
}


def wetted_area(W_TO, airplane_type: str):
    """Statistical wetted area (m^2) from take-off weight (N) (``WettedArea``)."""
    c, d = WETTED_AREA_COEFFS[airplane_type]
    return 10 ** (c + d * np.log10(np.asarray(W_TO) / LBF)) * FT**2


def equivalent_parasite_area(c_f, S_wet):
    """Equivalent parasite area ``f = c_f S_wet`` (``EquivalentParasiteArea``).

    The MATLAB fitted a polynomial through Roskam's (c_f, a) table; that
    table is exactly ``a = log10(c_f)`` with ``b = 1``.
    """
    return c_f * np.asarray(S_wet)


def zero_lift_drag(f, S):
    """C_D0 = f / S (``DragCoefficientZeroLift``)."""
    return f / S


def wetted_area_planform(S_exposed, tc_root, tc_tip, c_root, c_tip):
    """Class I wetted area of a lifting surface (``WettedArea_Planform_Class_I``).

    ``2 S_exp (1 + 0.25 (t/c)_r (1 + tau lambda) / (1 + lambda))`` with
    ``tau = (t/c)_t/(t/c)_r``. The MATLAB grouping of these terms was garbled
    (``4 * 1 + c_t/c_r``, inverted tau); this is the textbook form.
    """
    tau = tc_tip / tc_root
    lam = c_tip / c_root
    return 2 * S_exposed * (1 + 0.25 * tc_root * (1 + tau * lam) / (1 + lam))


def wetted_area_fuselage(D_f, l_f, cylindrical_mid_section: bool, l_nose=None):
    """Class I fuselage wetted area (``WettedArea_Fuselage_Class_I``).

    Textbook forms: ``pi D l (1 - 2/f)^(2/3) (1 + 1/f^2)`` with a cylindrical
    mid-section, ``pi D l (0.5 + 0.135 l_n/l)^(2/3) (1.015 + 0.3/f^1.5)``
    streamlined (f = l/D). The MATLAB misplaced parentheses in the first and
    had 1.35 for 0.135 in the second.
    """
    fineness = l_f / D_f
    if cylindrical_mid_section:
        return np.pi * D_f * l_f * (1 - 2 / fineness) ** (2 / 3) * (1 + 1 / fineness**2)
    return np.pi * D_f * l_f * (0.5 + 0.135 * l_nose / l_f) ** (2 / 3) * (1.015 + 0.3 / fineness**1.5)


# Engines


def piston_power_ratio(rho_ratio=1.0, T_ratio=1.0, mass_ratio=1.0, piston_area_ratio=1.0, speed_ratio=1.0, energy_ratio=1.0):
    """Barrett's IC-engine predictor: BHP_1/BHP_0 as a product of ratios
    (``Internal_Combustion_Engine_Predictor_Equation``).

    ``rho_ratio = rho1/rho0``, ``T_ratio = T1/T0`` (enters as sqrt),
    ``mass_ratio`` is ``(1/29 + f/mdot_f + h/18)_0 / (...)_1``, ``energy_ratio`` is
    ``(f Q_c eta_th eta_me)_1 / (...)_0``.
    """
    return rho_ratio * mass_ratio * piston_area_ratio * speed_ratio * np.sqrt(T_ratio) * energy_ratio


def piston_power_lapse(P_rated, h_rated, h):
    """Available shaft power at altitude ``h`` for an engine flat-rated at ``h_rated`` (``Power_Lapse_With_Alt_IC_engine``)."""
    a, b = isa(h_rated), isa(h)
    return P_rated * piston_power_ratio(rho_ratio=b.rho / a.rho, T_ratio=b.T / a.T)


# Climb


def climb_cl_max_rate(C_D0, A, e):
    """C_L for best rate of climb / minimum sink, ``sqrt(3 C_D0 pi A e)`` (``LiftCoefficent_RateOfClimb_Max``)."""
    return np.sqrt(3 * C_D0 * np.pi * A * e)


def climb_cd_max_rate(C_D0):
    """C_D at best rate of climb, ``4 C_D0`` (``DragCoefficent_RateOfClimb_Max``)."""
    return 4 * C_D0


def climb_index_max(A, e, C_D0):
    """Maximum ``C_L^1.5 / C_D = 1.345 (A e)^0.75 / C_D0^0.25`` (``ClimbIndex``/``CL32_CD_max``)."""
    return 1.345 * (A * e) ** 0.75 / C_D0**0.25


def power_for_climb(RC, eta_p, W, wing_loading, sigma, C_D0, A, e):
    """Shaft power (W) to climb at ``RC`` (m/s), Roskam Part I eq. 3.37 (``P_Climb``/``W2P_Climb_Propeller``).

    ``P[hp] = W/eta_p (RC[fpm]/33000 + sqrt(W/S[psf]) / (19 climb_index sqrt(sigma)))``.
    """
    W_lbf = np.asarray(W) / LBF
    rc_fpm = np.asarray(RC) / FPM
    ws = np.asarray(wing_loading) / PSF
    hp = W_lbf / eta_p * (rc_fpm / 33000 + np.sqrt(ws) / (19 * climb_index_max(A, e, C_D0) * np.sqrt(sigma)))
    return hp * HP


def power_loading_climb(RC, eta_p, W, wing_loading, sigma, C_D0, A, e):
    """Climb power loading W/P (N/W) (``W2P_Climb_Propeller``)."""
    return np.asarray(W) / power_for_climb(RC, eta_p, W, wing_loading, sigma, C_D0, A, e)


def power_required(rho, V, W, S, C_D0, A, e, n=1.0):
    """Aerodynamic power to fly at ``V`` with load factor ``n`` (``Required_Power_Subsonic``).

    The MATLAB squared the weight in ``C_L``; it should be ``C_L = 2 n W / (rho V^2 S)``.
    """
    V = np.asarray(V, dtype=float)
    C_L = 2 * n * W / (rho * V**2 * S)
    return 0.5 * rho * V**3 * S * parabolic_drag_polar(C_L, C_D0, A, e)


def electric_power(W, V, L_D, eta_p, eta_elec):
    """Battery power for a flight segment ``W V / (L/D eta_p eta_elec)`` (``powerest``)."""
    return W * np.asarray(V) / (L_D * eta_p * eta_elec)


def rate_of_climb_for_time(t_climb, h, h_abs):
    """Sea-level rate of climb to reach ``h`` in ``t_climb`` with absolute ceiling ``h_abs`` (Roskam eq. 3.33)."""
    return h_abs / t_climb * np.log(1 / (1 - h / h_abs))


def _steep_climb_sin_gamma(T2W, L_D_max):
    P = L_D_max**2 / (1 + L_D_max**2)
    return T2W * (P - np.sqrt(P**2 - P + 1 / ((1 + L_D_max**2) * T2W**2)))


# Constraint-diagram lines (jets, FAR 25) - all return T/W or W/S


def t2w_climb_gradient(CGR, L_D, n_engines, one_engine_inoperative):
    """T/W for a climb gradient, scaled up for one engine out (``T2W_ClimbGradient_FAR25``)."""
    oei = np.asarray(one_engine_inoperative, dtype=bool)
    factor = np.where(oei, n_engines / (n_engines - 1), 1.0)
    return (1 / np.asarray(L_D) + CGR) * factor


def t2w_takeoff_far25(field_length, wing_loading_TO, C_L_max_TO, sigma=1.0):
    """T/W for a FAR 25 take-off field length (m) (``T2W_Takeoff_FAR25``): ``37.5 (W/S) / (sigma C_Lmax S_TOFL)``."""
    return 37.5 * (np.asarray(wing_loading_TO) / PSF) / (sigma * C_L_max_TO * field_length / FT)


def t2w_cruise_speed(V, A, e, C_D0, wing_loading_TO, weight_ratio=1.0, sigma=1.0):
    """T/W to cruise at ``V`` (``T2W_CruiseSpeed``)."""
    q = 0.5 * RHO0 * sigma * V**2
    ws = np.asarray(wing_loading_TO) * weight_ratio
    return C_D0 * q / ws + ws / (q * np.pi * A * e)


def t2w_time_to_climb(t_climb, h, h_abs, A, e, C_D0, wing_loading, sigma=1.0, steep: bool = False):
    """T/W to meet a time-to-climb requirement (``T2W_TimeToClimb``). NaN if unattainable."""
    rho = RHO0 * sigma
    L_D = max_lift_to_drag(C_D0, A, e)
    V = np.sqrt(2 * np.asarray(wing_loading) / (rho * np.sqrt(C_D0 * np.pi * A * e)))
    RC = rate_of_climb_for_time(t_climb, h, h_abs)
    if not steep:
        return RC * L_D / V + 1

    tw_max = np.sqrt(1 + L_D**2) / L_D  # beyond this the steep-climb formula has no real root

    def solve(v):
        tw = np.linspace(1e-6, tw_max, 4001)
        f = _steep_climb_sin_gamma(tw, L_D) * v - RC
        k = np.nonzero(np.diff(np.sign(f)))[0]
        if not len(k):
            return np.nan
        return brentq(lambda t: _steep_climb_sin_gamma(t, L_D) * v - RC, tw[k[0]], tw[k[0] + 1])

    return np.vectorize(solve)(V)


def w2s_time_to_climb(t_climb, h, h_abs, A, e, C_D0, T2W, sigma=1.0, steep: bool = False):
    """W/S that meets a time-to-climb requirement at a given T/W (``W2S_TimeToClimb``)."""
    rho = RHO0 * sigma
    L_D = max_lift_to_drag(C_D0, A, e)
    RC = rate_of_climb_for_time(t_climb, h, h_abs)
    climb_ratio = _steep_climb_sin_gamma(T2W, L_D) if steep else (np.asarray(T2W) - 1) / L_D
    return (RC / climb_ratio) ** 2 * rho * np.sqrt(C_D0 * np.pi * A * e) / 2


def w2s_landing_far25(field_length, landing_weight_ratio, C_L_max_L, sigma=1.0):
    """Maximum take-off W/S for a FAR 25 landing field length (m) (``W2S_Landing_FAR25``).

    ``S_FL[ft] = 0.3 V_A[kt]^2``, ``V_A = 1.3 V_SL``.
    """
    from ..units import KT

    V_A = np.sqrt(np.asarray(field_length) / FT / 0.3) * KT
    V_S = V_A / 1.3
    return 0.5 * RHO0 * sigma * V_S**2 * C_L_max_L / landing_weight_ratio


# FAR 25 climb-gradient requirements by engine count:
# (initial, transition, second segment, en-route, balked landing AEO, balked landing OEI)
FAR25_CLIMB_GRADIENTS = {
    2: (0.012, 0.0, 0.024, 0.012, 0.032, 0.021),
    3: (0.015, 0.003, 0.027, 0.015, 0.032, 0.024),
    4: (0.017, 0.005, 0.030, 0.017, 0.032, 0.027),
}


def t2w_climb_far25(
    n_engines,
    wing_loading_clean,
    wing_loading_TO,
    wing_loading_L,
    C_L_max_clean,
    C_L_max_TO,
    C_L_max_L,
    drag_polars,
    W_MGTO,
    W_L_max,
    hot_day_thrust_ratio=0.8,
):
    """T/W (referred to MGTO, sea-level standard) for the six FAR 25 climb cases (``T2W_Climb_FAR25``).

    ``drag_polars`` is a sequence of six callables ``C_D(C_L)`` for:
    clean, clean+gear, take-off flaps, take-off flaps+gear, landing flaps,
    landing flaps+gear. Returns the six T/W values in the order of
    ``FAR25_CLIMB_GRADIENTS``.
    """
    CGR = np.array(FAR25_CLIMB_GRADIENTS[n_engines])
    C_L = np.array([C_L_max_TO, C_L_max_TO, C_L_max_TO, C_L_max_clean, C_L_max_L, C_L_max_L])
    polar_index = [2, 3, 2, 0, 5, 5]
    C_D = np.array([drag_polars[i](cl) for i, cl in zip(polar_index, C_L)])
    oei = np.array([1, 1, 1, 1, 0, 1], dtype=bool)
    W = np.array([W_MGTO] * 4 + [W_L_max] * 2)
    T2W = t2w_climb_gradient(CGR, C_L / C_D, n_engines, oei) / hot_day_thrust_ratio
    return T2W * W / W_MGTO


# Gliders


def dive_speed(wing_loading, C_D0, A, e, rho, gamma):
    """Steady 1 g dive speed at flight-path angle ``gamma`` (rad) (``SpeedDive``)."""
    return np.sqrt(
        2 * wing_loading * (np.sin(gamma) - 2 / rho * wing_loading * C_D0 * np.cos(gamma) ** 2 / (np.pi * A * e))
        / (rho * C_D0)
    )


def airbrake_drag_coefficient(area_ratio):
    """Schempp-Hirth airbrake increment ``0.3 S_AB/S`` (Panjo, Sailplane Design) (``DragCoefficient_Airbrakes``)."""
    return 0.3 * np.asarray(area_ratio)


def w2s_for_lift_to_drag(L_D, V, C_D0, C_L_max, A, e, rho):
    """Wing-loading band ``(low, high)`` that gives at least ``L_D`` at speed ``V`` (``W2S_LiftToDragAtSpeed``).

    The upper bound is capped at the stall wing loading ``C_Lmax rho V^2 / 2``.
    NaN where the target exceeds (L/D)_max.
    """
    with np.errstate(invalid="ignore"):
        root = np.sqrt((np.pi * A * e - 4 * C_D0 * L_D**2) / (A * e))
    k = A * V**2 * e * rho * np.sqrt(np.pi) / (4 * L_D)
    low = k * (np.sqrt(np.pi) - root)
    high = np.minimum(k * (np.sqrt(np.pi) + root), C_L_max * rho * V**2 / 2)
    return low, high


def battery_mass(power, energy, specific_power, specific_energy_wh):
    """Battery mass (kg) sized by whichever of power (W/kg) or energy (Wh/kg) governs (AE 722 weight sizing)."""
    return np.maximum(power / specific_power, energy / 3600 / specific_energy_wh)


def climb_power_far(RC, eta_p, W, S, rho, C_D0, A, e):
    """Climb power from ``RC W / eta + 0.8776 sqrt(W/S / (rho C_D0)) W / (eta C_L^1.5/C_D)`` (AE 722 ``Climb_Power``)."""
    C_L = climb_cl_max_rate(C_D0, A, e)
    C_D = climb_cd_max_rate(C_D0)
    return RC * W / eta_p + 0.8776 * np.sqrt(W / S / (rho * C_D0)) * W / (eta_p * C_L**1.5 / C_D)


def slug_density(rho):
    """kg/m^3 -> slug/ft^3."""
    return rho / SLUG_PER_FT3


# FAR 23 / CS-22 propeller-aircraft constraints (AE 722 Archytas constraint diagram)


def landing_wing_loading_far23(landing_distance, C_L_max_L, landing_weight_ratio=1.0, rho=RHO0):
    """Max take-off W/S (N/m^2) for a FAR 23 landing distance (m) (``Landing_Requirements_Reida``).

    Ground run ``S_LG = S_L / 1.938``; ``S_LG[ft] = 0.265 V_SL[kt]^2``.
    """
    from ..units import KT

    V_SL = np.sqrt(np.asarray(landing_distance) / 1.938 / FT / 0.265) * KT
    return 0.5 * rho * V_SL**2 * np.asarray(C_L_max_L) / landing_weight_ratio


def takeoff_power_loading_far23(TOP23, wing_loading, C_L_max_TO, sigma=1.0):
    """W/P (N/W) from Roskam's FAR 23 take-off parameter (lbf^2/(ft^2 hp)) (``Takeoff_Requirements_Reworked``).

    ``(W/P)[lbf/hp] = TOP23 sigma C_LmaxTO / (W/S)[psf]``.
    """
    wp = TOP23 * sigma * np.asarray(C_L_max_TO) / (np.asarray(wing_loading) / PSF)
    return wp * LBF / HP


def sea_level_rate_of_climb(RC_at_ceiling, h, h_abs):
    """Sea-level rate of climb implied by ``RC_at_ceiling`` at ``h`` with a linear lapse to ``h_abs``."""
    return RC_at_ceiling / (1 - h / h_abs)
