"""Gas-turbine cycle analysis (AE 573): on-design parametric analysis with component losses.

Replaces the equation lists of the AE 573 scripts (``HomeWork``, ``Quiz_2``,
``HW_7``, ``HW_8_problem_2``, ``H_8_part_3``, ``HW_10``, ``EXAM_1_simple``,
``Exam_2_part_1``, ``Exam_2_part_3``, ``Pratice_quiz_turbofan``,
``FinalExam_1``, ``FinalExam_3``), which mixed up ratios (``tau_c_lp = Tt_2 /
Tt_2_5``, ``pi_d = pt_0 / pt_2``), dropped parentheses (``cp_4*Tt_4 /
cp_0*T_0``, ``(p_0/pt_45)^(gamma-1)/gamma``), used ``M`` for ``M^2`` and
held relations that contradicted each other, so the relation solver picked
whichever it reached first. Here every engine is written out explicitly
in the order of Mattingly, *Elements of Propulsion*: free stream, diffuser,
compressor/fan (polytropic efficiencies), burner (``f`` from the energy
balance with burner efficiency), turbine (from the shaft power balance
with mechanical efficiency), optional mixer and afterburner, nozzles.

Cold and hot sections can have different ``gamma`` and ``c_p``
(:class:`Gas`); ``R = c_p (gamma - 1) / gamma`` for each. Specific thrust
includes the pressure term ``R T9 / V9 (1 - p0 / p9)`` when a nozzle is not
fully expanded. Station names follow Mattingly: 0 free stream, 2 fan/
compressor face, 13 fan exit, 2.5 (``"25"``) LPC exit, 3 compressor exit, 4
burner exit, 4.5 (``"45"``) HPT exit, 5 turbine exit, 6A mixer exit, 7
afterburner exit, 9 and 19 core and fan nozzle exits.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class Gas:
    gamma: float = 1.4
    cp: float = 1004.0

    @property
    def R(self) -> float:
        return self.cp * (self.gamma - 1) / self.gamma


AIR = Gas(1.4, 1004.0)
COMBUSTION_GAS = Gas(1.33, 1156.0)


@dataclass
class Station:
    pt: float
    Tt: float


@dataclass
class NozzleExit:
    p: float
    T: float
    M: float
    V: float
    gas: Gas
    choked: bool = False

    def pressure_thrust_per_mass(self, p0) -> float:
        """``(p9 - p0) A9 / mdot9 = R T9 / V9 (1 - p0 / p9)``."""
        return self.gas.R * self.T / self.V * (1 - p0 / self.p) if self.V > 0 else 0.0

    def effective_velocity(self, p0) -> float:
        """Exhaust velocity including the pressure thrust, ``V9 + R T9 / V9 (1 - p0 / p9)``."""
        return self.V + self.pressure_thrust_per_mass(p0)


@dataclass
class CycleResult:
    """Performance per unit mass flow. ``specific_thrust`` is ``F / mdot_0`` with
    ``mdot_0`` the total (core + bypass) air flow; ``tsfc`` is in kg/(N s)."""

    specific_thrust: float
    tsfc: float
    f: float
    V0: float
    stations: dict
    exits: dict
    f_ab: float = 0.0
    alpha: float = 0.0
    eta_thermal: float = np.nan
    eta_propulsive: float = np.nan
    eta_overall: float = np.nan
    extra: dict = field(default_factory=dict)

    def thrust(self, mdot0) -> float:
        return self.specific_thrust * mdot0

    def summary(self) -> str:
        lines = [f"F/mdot0 = {self.specific_thrust:.2f} N s/kg, TSFC = {self.tsfc * 1e6:.3f} mg/(N s)",
                 f"f = {self.f:.5f}" + (f", f_AB = {self.f_ab:.5f}" if self.f_ab else "")
                 + (f", alpha = {self.alpha:.4f}" if self.alpha else ""),
                 f"eta_th = {self.eta_thermal:.4f}, eta_p = {self.eta_propulsive:.4f}, eta_o = {self.eta_overall:.4f}"]
        for k, s in self.stations.items():
            lines.append(f"  station {k:>3s}: pt = {s.pt:12.1f} Pa, Tt = {s.Tt:8.2f} K")
        for k, e in self.exits.items():
            lines.append(f"  exit {k:>3s}: p = {e.p:10.1f} Pa, T = {e.T:8.2f} K, M = {e.M:.4f}, V = {e.V:8.2f} m/s"
                         + (" (choked)" if e.choked else ""))
        for k, v in self.extra.items():
            lines.append(f"  {k} = {v:.6g}")
        return "\n".join(lines)


# --- components ---------------------------------------------------------------------------------

def total_temperature_ratio(M, gas: Gas = AIR):
    return 1 + (gas.gamma - 1) / 2 * np.asarray(M) ** 2


def total_pressure_ratio(M, gas: Gas = AIR):
    return total_temperature_ratio(M, gas) ** (gas.gamma / (gas.gamma - 1))


def freestream(M0, T0, p0, gas: Gas = AIR):
    """``(V0, Station0, tau_r, pi_r)``."""
    tau_r = float(total_temperature_ratio(M0, gas))
    pi_r = tau_r ** (gas.gamma / (gas.gamma - 1))
    V0 = M0 * np.sqrt(gas.gamma * gas.R * T0)
    return V0, Station(p0 * pi_r, T0 * tau_r), tau_r, pi_r


def compressor(pi, e=1.0, gas: Gas = AIR):
    """``(tau, eta)`` of a compressor or fan with polytropic efficiency ``e``."""
    g = gas.gamma
    tau = pi ** ((g - 1) / (g * e))
    return tau, (pi ** ((g - 1) / g) - 1) / (tau - 1)


def turbine(tau, e=1.0, gas: Gas = COMBUSTION_GAS):
    """``(pi, eta)`` of a turbine with temperature ratio ``tau`` and polytropic efficiency ``e``."""
    g = gas.gamma
    pi = tau ** (g / ((g - 1) * e))
    return pi, (1 - tau) / (1 - pi ** ((g - 1) / g))


def compressor_efficiencies(pi, tau, gas: Gas = AIR):
    """Measured ``pi`` and ``tau`` -> ``(eta_adiabatic, e_polytropic)``."""
    g = gas.gamma
    return (pi ** ((g - 1) / g) - 1) / (tau - 1), (g - 1) / g * np.log(pi) / np.log(tau)


def turbine_efficiencies(pi, tau, gas: Gas = COMBUSTION_GAS):
    """Measured ``pi`` and ``tau`` (both < 1) -> ``(eta_adiabatic, e_polytropic)``."""
    g = gas.gamma
    return (1 - tau) / (1 - pi ** ((g - 1) / g)), g / (g - 1) * np.log(tau) / np.log(pi)


def burner_fuel_air_ratio(Tt_in, Tt_out, Qr, eta_b=1.0, gas_in: Gas = AIR, gas_out: Gas = COMBUSTION_GAS):
    """``f`` from ``(1 + f) cp_out Tt_out - cp_in Tt_in = eta_b f Qr``."""
    return (gas_out.cp * Tt_out - gas_in.cp * Tt_in) / (eta_b * Qr - gas_out.cp * Tt_out)


def nozzle(pt, Tt, p0, gas: Gas, convergent=False, p_exit=None) -> NozzleExit:
    """Exit state of a nozzle fed at ``pt``, ``Tt`` (losses already in ``pt``).

    By default the nozzle is fully expanded to ``p_exit`` (``p0`` if not
    given). ``convergent=True`` gives a convergent nozzle: exit at ``p0`` if
    the flow is subsonic, otherwise choked at ``M = 1`` with ``p > p0``.
    """
    g = gas.gamma
    if convergent:
        p_star = pt / ((g + 1) / 2) ** (g / (g - 1))
        p, choked = (p_star, True) if p_star > p0 else (p0, False)
    else:
        p, choked = (p0 if p_exit is None else p_exit), False
    ratio = max((pt / p) ** ((g - 1) / g), 1.0)
    M = np.sqrt(2 / (g - 1) * (ratio - 1))
    T = Tt / ratio
    return NozzleExit(p, T, M, M * np.sqrt(g * gas.R * T), gas, choked)


def _efficiencies(result: CycleResult, kinetic_power, fuel_power):
    """Thermal, propulsive and overall efficiency from the kinetic energy rise per ``mdot_0``.

    The kinetic energy uses the effective exhaust velocities (pressure thrust
    included), so that ``eta_p = 2 V0 / (V_eff + V0)`` stays below 1 for
    nozzles that are not fully expanded.
    """
    result.eta_thermal = kinetic_power / fuel_power
    result.eta_propulsive = result.specific_thrust * result.V0 / kinetic_power if result.V0 > 0 else 0.0
    result.eta_overall = result.eta_thermal * result.eta_propulsive
    return result


# --- engines ------------------------------------------------------------------------------------

def turbojet(M0, T0, p0, pi_c, Tt4, Qr, *, gas_c: Gas = AIR, gas_t: Gas = COMBUSTION_GAS, pi_d=1.0, pi_b=1.0,
             pi_n=1.0, e_c=1.0, e_t=1.0, eta_b=1.0, eta_m=1.0, Tt7=None, gas_ab: Gas | None = None, pi_ab=1.0,
             eta_ab=1.0, convergent=False, p9=None) -> CycleResult:
    """Single-spool turbojet, with an afterburner when ``Tt7`` is given (Mattingly 7.2)."""
    V0, st0, tau_r, pi_r = freestream(M0, T0, p0, gas_c)
    st = {"0": st0, "2": Station(st0.pt * pi_d, st0.Tt)}
    tau_c, eta_c = compressor(pi_c, e_c, gas_c)
    st["3"] = Station(st["2"].pt * pi_c, st["2"].Tt * tau_c)
    f = burner_fuel_air_ratio(st["3"].Tt, Tt4, Qr, eta_b, gas_c, gas_t)
    st["4"] = Station(st["3"].pt * pi_b, Tt4)
    Tt5 = Tt4 - gas_c.cp * (st["3"].Tt - st["2"].Tt) / (eta_m * (1 + f) * gas_t.cp)
    tau_t = Tt5 / Tt4
    pi_t, eta_t = turbine(tau_t, e_t, gas_t)
    st["5"] = Station(st["4"].pt * pi_t, Tt5)
    f_ab, gas_n, last = 0.0, gas_t, st["5"]
    if Tt7 is not None:
        gas_n = gas_ab or gas_t
        f_ab = (1 + f) * (gas_n.cp * Tt7 - gas_t.cp * Tt5) / (eta_ab * Qr - gas_n.cp * Tt7)
        last = st["7"] = Station(st["5"].pt * pi_ab, Tt7)
    st["9"] = Station(last.pt * pi_n, last.Tt)
    ex = nozzle(st["9"].pt, st["9"].Tt, p0, gas_n, convergent, p9)
    m9 = 1 + f + f_ab
    F = m9 * ex.V - V0 + m9 * ex.pressure_thrust_per_mass(p0)
    res = CycleResult(F, (f + f_ab) / F, f, V0, st, {"9": ex}, f_ab=f_ab,
                      extra={"tau_r": tau_r, "pi_r": pi_r, "tau_c": tau_c, "eta_c": eta_c, "tau_t": tau_t,
                             "pi_t": pi_t, "eta_t": eta_t})
    return _efficiencies(res, (m9 * ex.effective_velocity(p0) ** 2 - V0**2) / 2, (f + f_ab) * Qr)


def turbofan(M0, T0, p0, pi_c, pi_f, alpha, Tt4, Qr, *, gas_c: Gas = AIR, gas_t: Gas = COMBUSTION_GAS, pi_d=1.0,
             pi_b=1.0, pi_n=1.0, pi_fn=1.0, e_c=1.0, e_f=1.0, e_t=1.0, eta_b=1.0, eta_m=1.0, convergent=False,
             convergent_fan=None) -> CycleResult:
    """Separate-exhaust turbofan (Mattingly 7.3). ``pi_c`` is the overall core
    pressure ratio (fan root included); one turbine drives fan and compressor:
    ``eta_m (1 + f) cp_t (Tt4 - Tt5) = cp_c [(Tt3 - Tt2) + alpha (Tt13 - Tt2)]``.
    ``specific_thrust`` is per total air flow ``(1 + alpha) mdot_core``.
    """
    V0, st0, tau_r, pi_r = freestream(M0, T0, p0, gas_c)
    st = {"0": st0, "2": Station(st0.pt * pi_d, st0.Tt)}
    tau_f, eta_f = compressor(pi_f, e_f, gas_c)
    tau_c, eta_c = compressor(pi_c, e_c, gas_c)
    st["13"] = Station(st["2"].pt * pi_f, st["2"].Tt * tau_f)
    st["3"] = Station(st["2"].pt * pi_c, st["2"].Tt * tau_c)
    f = burner_fuel_air_ratio(st["3"].Tt, Tt4, Qr, eta_b, gas_c, gas_t)
    st["4"] = Station(st["3"].pt * pi_b, Tt4)
    work = gas_c.cp * ((st["3"].Tt - st["2"].Tt) + alpha * (st["13"].Tt - st["2"].Tt))
    Tt5 = Tt4 - work / (eta_m * (1 + f) * gas_t.cp)
    if Tt5 <= 0:
        raise ValueError("the turbine cannot supply the fan and compressor work (Tt5 <= 0)")
    tau_t = Tt5 / Tt4
    pi_t, eta_t = turbine(tau_t, e_t, gas_t)
    st["5"] = Station(st["4"].pt * pi_t, Tt5)
    st["9"] = Station(st["5"].pt * pi_n, Tt5)
    st["19"] = Station(st["13"].pt * pi_fn, st["13"].Tt)
    ex9 = nozzle(st["9"].pt, Tt5, p0, gas_t, convergent)
    ex19 = nozzle(st["19"].pt, st["19"].Tt, p0, gas_c, convergent if convergent_fan is None else convergent_fan)
    F_core = (1 + f) * ex9.V - V0 + (1 + f) * ex9.pressure_thrust_per_mass(p0)
    F_fan = ex19.V - V0 + ex19.pressure_thrust_per_mass(p0)
    F = (F_core + alpha * F_fan) / (1 + alpha)
    res = CycleResult(F, f / ((1 + alpha) * F), f, V0, st, {"9": ex9, "19": ex19}, alpha=alpha,
                      extra={"F_core/mdot_core": F_core, "F_fan/mdot_fan": F_fan, "tau_r": tau_r, "pi_r": pi_r,
                             "tau_f": tau_f, "eta_f": eta_f, "tau_c": tau_c, "eta_c": eta_c, "tau_t": tau_t,
                             "pi_t": pi_t, "eta_t": eta_t})
    ke = ((1 + f) * ex9.effective_velocity(p0) ** 2 + alpha * ex19.effective_velocity(p0) ** 2
          - (1 + alpha) * V0**2) / 2 / (1 + alpha)
    return _efficiencies(res, ke, f * Qr / (1 + alpha))


def _mfp(M, gas: Gas):
    """Mass-flow parameter ``mdot sqrt(Tt) / (pt A)``."""
    g = gas.gamma
    return M * np.sqrt(g / gas.R) * total_temperature_ratio(M, gas) ** (-(g + 1) / (2 * (g - 1)))


def constant_area_mixer(m_core, core: Station, gas_core: Gas, M_core, m_bypass, bypass: Station, gas_bypass: Gas):
    """Ideal constant-area mixer (Mattingly 7.4): equal static pressure at entry,
    conservation of mass, energy and impulse ``p A (1 + gamma M^2)``.

    Returns ``(Station6A, Gas6A, M6A, M_bypass)``; the exit total pressure
    over the core total pressure is the ideal mixer ratio ``pi_M``.
    """
    g6 = gas_core.gamma
    p6 = core.pt / total_pressure_ratio(M_core, gas_core)
    g16 = gas_bypass.gamma
    ratio16 = (bypass.pt / p6) ** ((g16 - 1) / g16)
    if ratio16 < 1:
        raise ValueError("bypass total pressure below the core static pressure: the streams cannot mix")
    M16 = np.sqrt(2 / (g16 - 1) * (ratio16 - 1))
    A6 = m_core * np.sqrt(core.Tt) / (core.pt * _mfp(M_core, gas_core))
    A16 = m_bypass * np.sqrt(bypass.Tt) / (bypass.pt * _mfp(M16, gas_bypass))
    m6A = m_core + m_bypass
    cp6A = (m_core * gas_core.cp + m_bypass * gas_bypass.cp) / m6A
    R6A = (m_core * gas_core.R + m_bypass * gas_bypass.R) / m6A
    gas6A = Gas(cp6A / (cp6A - R6A), cp6A)
    Tt6A = (m_core * gas_core.cp * core.Tt + m_bypass * gas_bypass.cp * bypass.Tt) / (m6A * cp6A)
    impulse = p6 * A6 * (1 + g6 * M_core**2) + p6 * A16 * (1 + g16 * M16**2)
    A6A = A6 + A16
    g = gas6A.gamma
    target = m6A**2 * gas6A.R * Tt6A / (g * impulse**2)

    def phi(M):
        return M**2 * (1 + (g - 1) / 2 * M**2) / (1 + g * M**2) ** 2 - target

    M6A = brentq(phi, 1e-6, 1.0)  # subsonic root
    p6A = impulse / (A6A * (1 + g * M6A**2))
    return Station(p6A * total_pressure_ratio(M6A, gas6A), Tt6A), gas6A, M6A, M16


def mixed_turbofan(M0, T0, p0, pi_c, pi_f, Tt4, Qr, M6, *, gas_c: Gas = AIR, gas_t: Gas = COMBUSTION_GAS,
                   pi_d=1.0, pi_b=1.0, pi_fd=1.0, pi_m_friction=1.0, pi_n=1.0, e_c=1.0, e_f=1.0, e_t=1.0,
                   eta_b=1.0, eta_m=1.0, Tt7=None, gas_ab: Gas | None = None, pi_ab=1.0, eta_ab=1.0,
                   pt_ratio_16_6=1.0, convergent=False, p9=None) -> CycleResult:
    """Mixed-flow turbofan with optional afterburner (Mattingly 7.4).

    The bypass ratio follows from matching the mixer inlet total pressures,
    ``pt16 = pt_ratio_16_6 * pt6`` with ``pt16 = pt2 pi_f pi_fd`` and ``pt6 =
    pt2 pi_c pi_b pi_t``: that fixes ``pi_t`` (``pi_c`` is the overall ratio),
    hence ``Tt5``, and the turbine power balance then gives ``alpha``. The
    mixer is the ideal constant-area mixer at core Mach ``M6`` times
    ``pi_m_friction``.
    """
    V0, st0, tau_r, pi_r = freestream(M0, T0, p0, gas_c)
    st = {"0": st0, "2": Station(st0.pt * pi_d, st0.Tt)}
    tau_f, eta_f = compressor(pi_f, e_f, gas_c)
    tau_c, eta_c = compressor(pi_c, e_c, gas_c)
    st["13"] = Station(st["2"].pt * pi_f, st["2"].Tt * tau_f)
    st["3"] = Station(st["2"].pt * pi_c, st["2"].Tt * tau_c)
    f = burner_fuel_air_ratio(st["3"].Tt, Tt4, Qr, eta_b, gas_c, gas_t)
    st["4"] = Station(st["3"].pt * pi_b, Tt4)
    pi_t = pi_f * pi_fd / (pi_c * pi_b * pt_ratio_16_6)
    tau_t = pi_t ** ((gas_t.gamma - 1) * e_t / gas_t.gamma)
    _, eta_t = turbine(tau_t, e_t, gas_t)
    st["5"] = Station(st["4"].pt * pi_t, Tt4 * tau_t)
    # eta_m (1 + f) cp_t (Tt4 - Tt5) = cp_c [(Tt3 - Tt2) + alpha (Tt13 - Tt2)]
    alpha = ((eta_m * (1 + f) * gas_t.cp * (Tt4 - st["5"].Tt) - gas_c.cp * (st["3"].Tt - st["2"].Tt))
             / (gas_c.cp * (st["13"].Tt - st["2"].Tt)))
    if alpha <= 0:
        raise ValueError(f"no positive bypass ratio matches the mixer pressures (alpha = {alpha:.4g})")
    st["16"] = Station(st["13"].pt * pi_fd, st["13"].Tt)
    st6A, gas6A, M6A, M16 = constant_area_mixer(1 + f, st["5"], gas_t, M6, alpha, st["16"], gas_c)
    pi_m_ideal = st6A.pt / st["5"].pt
    st6A.pt *= pi_m_friction
    st["6A"] = st6A
    f_ab, gas_n, last = 0.0, gas6A, st6A
    if Tt7 is not None:
        gas_n = gas_ab or gas6A
        # per unit core air: (1 + f + alpha) and the afterburner fuel f_ab
        f_ab = (1 + f + alpha) * (gas_n.cp * Tt7 - gas6A.cp * st6A.Tt) / (eta_ab * Qr - gas_n.cp * Tt7)
        last = st["7"] = Station(st6A.pt * pi_ab, Tt7)
    st["9"] = Station(last.pt * pi_n, last.Tt)
    ex = nozzle(st["9"].pt, st["9"].Tt, p0, gas_n, convergent, p9)
    m9 = 1 + f + alpha + f_ab  # per unit core air
    F = (m9 * ex.V - (1 + alpha) * V0 + m9 * ex.pressure_thrust_per_mass(p0)) / (1 + alpha)
    res = CycleResult(F, (f + f_ab) / ((1 + alpha) * F), f, V0, st, {"9": ex}, f_ab=f_ab, alpha=alpha,
                      extra={"tau_r": tau_r, "pi_r": pi_r, "tau_f": tau_f, "tau_c": tau_c, "tau_t": tau_t,
                             "pi_t": pi_t, "eta_t": eta_t, "pi_M_ideal": pi_m_ideal, "M6A": M6A, "M16": M16,
                             "gamma_6A": gas6A.gamma, "cp_6A": gas6A.cp})
    ke = (m9 * ex.effective_velocity(p0) ** 2 - (1 + alpha) * V0**2) / 2 / (1 + alpha)
    return _efficiencies(res, ke, (f + f_ab) * Qr / (1 + alpha))


def turboprop(M0, T0, p0, pi_c, Tt4, Qr, work_split, *, gas_c: Gas = AIR, gas_t: Gas = COMBUSTION_GAS, pi_d=1.0,
              pi_b=1.0, e_c=1.0, e_hpt=1.0, eta_b=1.0, eta_m_hpt=1.0, eta_lpt=1.0, eta_m_lpt=1.0, eta_gb=1.0,
              eta_prop=1.0, eta_n=1.0) -> CycleResult:
    """Turboprop with the available expansion work split between the free turbine and the nozzle (Hill & Peterson).

    The HPT drives the compressor. From station 4.5 the isentropic work
    available down to ``p0`` is ``h_avail = cp_t Tt45 (1 - (p0 / pt45)^((g - 1)/g))``;
    the low-pressure turbine takes ``work_split`` of it with efficiency
    ``eta_lpt`` and drives the propeller through the gearbox, the nozzle
    expands the rest: ``V9^2 = 2 (1 - work_split) eta_n h_avail``.
    ``specific_thrust`` adds the propeller thrust ``eta_prop P_prop / V0``.
    """
    V0, st0, tau_r, pi_r = freestream(M0, T0, p0, gas_c)
    st = {"0": st0, "2": Station(st0.pt * pi_d, st0.Tt)}
    tau_c, eta_c = compressor(pi_c, e_c, gas_c)
    st["3"] = Station(st["2"].pt * pi_c, st["2"].Tt * tau_c)
    f = burner_fuel_air_ratio(st["3"].Tt, Tt4, Qr, eta_b, gas_c, gas_t)
    st["4"] = Station(st["3"].pt * pi_b, Tt4)
    Tt45 = Tt4 - gas_c.cp * (st["3"].Tt - st["2"].Tt) / (eta_m_hpt * (1 + f) * gas_t.cp)
    pi_hpt, eta_hpt = turbine(Tt45 / Tt4, e_hpt, gas_t)
    st["45"] = Station(st["4"].pt * pi_hpt, Tt45)
    g = gas_t.gamma
    h_avail = gas_t.cp * Tt45 * (1 - (p0 / st["45"].pt) ** ((g - 1) / g))
    w_lpt = eta_lpt * work_split * h_avail
    P_prop = eta_m_lpt * eta_gb * (1 + f) * w_lpt  # shaft power to the propeller per unit air flow
    V9 = np.sqrt(2 * (1 - work_split) * eta_n * h_avail)
    T9 = Tt45 - w_lpt / gas_t.cp - V9**2 / (2 * gas_t.cp)
    ex = NozzleExit(p0, T9, V9 / np.sqrt(g * gas_t.R * T9), V9, gas_t)
    st["5"] = Station(np.nan, Tt45 - w_lpt / gas_t.cp)
    F_core = (1 + f) * V9 - V0
    F_prop = eta_prop * P_prop / V0 if V0 > 0 else np.nan
    F = F_core + F_prop
    res = CycleResult(F, f / F, f, V0, st, {"9": ex},
                      extra={"P_prop/mdot0": P_prop, "F_prop/mdot0": F_prop, "F_core/mdot0": F_core,
                             "h_avail": h_avail, "tau_c": tau_c, "eta_c": eta_c, "pi_hpt": pi_hpt,
                             "eta_hpt": eta_hpt, "power_specific_fuel_consumption": f / P_prop})
    # thermal efficiency on shaft power plus jet kinetic energy, as in the MATLAB
    useful = P_prop + ((1 + f) * V9**2 - V0**2) / 2
    res.eta_thermal = useful / (f * Qr)
    res.eta_propulsive = F * V0 / useful if V0 > 0 else 0.0
    res.eta_overall = res.eta_thermal * res.eta_propulsive
    return res


# --- compressor stage ---------------------------------------------------------------------------

@dataclass
class CompressorStage:
    U: float
    ct1: float
    ct2: float
    work: float
    Tt3: float
    pt3: float
    pressure_ratio: float
    efficiency: float
    diffusion_factor_rotor: float
    diffusion_factor_stator: float
    M1: float
    M1_relative: float
    M2: float
    M3: float
    beta1: float
    beta2: float
    alpha2: float


def compressor_stage(radius, omega, cz, reaction, Tt1, pt1, solidity, loss_rotor=0.0, loss_stator=0.0,
                     ct1=0.0, gas: Gas = AIR) -> CompressorStage:
    """Mean-line axial compressor stage with constant axial velocity (repeating stage, ``c3 = c1``).

    Degree of reaction ``R = 1 - (ct1 + ct2) / (2 U)`` gives the rotor exit
    swirl; Euler's work ``U (ct2 - ct1)``; total-pressure loss coefficients
    ``omega = (pt_in - pt_out) / (pt_in - p_in)`` in the relative frame for
    the rotor and the absolute frame for the stator; Lieblein diffusion
    factors ``D = 1 - w2/w1 + |dct| / (2 sigma w1)``. Angles in degrees from
    the axial direction.
    """
    g, cp, R = gas.gamma, gas.cp, gas.R
    U = radius * omega
    ct2 = 2 * U * (1 - reaction) - ct1
    work = U * (ct2 - ct1)
    c1, c2 = np.hypot(cz, ct1), np.hypot(cz, ct2)
    w1, w2 = np.hypot(cz, U - ct1), np.hypot(cz, U - ct2)
    T1 = Tt1 - c1**2 / (2 * cp)
    p1 = pt1 * (T1 / Tt1) ** (g / (g - 1))
    # rotor, relative frame (rothalpy conserved, Ttr constant at the mean line)
    Ttr1 = T1 + w1**2 / (2 * cp)
    ptr1 = p1 * (Ttr1 / T1) ** (g / (g - 1))
    ptr2 = ptr1 - loss_rotor * (ptr1 - p1)
    Tt2 = Tt1 + work / cp
    T2 = Tt2 - c2**2 / (2 * cp)
    Ttr2 = T2 + w2**2 / (2 * cp)
    p2 = ptr2 * (T2 / Ttr2) ** (g / (g - 1))
    pt2 = p2 * (Tt2 / T2) ** (g / (g - 1))
    # stator back to c3 = c1
    pt3 = pt2 - loss_stator * (pt2 - p2)
    T3 = Tt2 - c1**2 / (2 * cp)
    pr = pt3 / pt1
    eta = (pr ** ((g - 1) / g) - 1) / (Tt2 / Tt1 - 1)
    D_r = 1 - w2 / w1 + abs(ct2 - ct1) / (2 * solidity * w1)
    D_s = 1 - c1 / c2 + abs(ct2 - ct1) / (2 * solidity * c2)
    a = lambda T: np.sqrt(g * R * T)  # noqa: E731
    return CompressorStage(U, ct1, ct2, work, Tt2, pt3, pr, eta, D_r, D_s, c1 / a(T1), w1 / a(T1), c2 / a(T2),
                           c1 / a(T3), np.degrees(np.arctan2(U - ct1, cz)), np.degrees(np.arctan2(U - ct2, cz)),
                           np.degrees(np.arctan2(ct2, cz)))
