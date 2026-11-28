"""AE 573 homework 8 problem 2 (``HW_8_problem_2``): component performance from turbofan test data.

Static sea-level test (V0 = 0, p0 = 101 325 Pa) of a separate-exhaust
turbofan: pt2 = 14.7, pt13 = 26, pt2.5 = 63, pt3 = 200, pt4 = 190,
pt5 = 28 psia; Tt2 = 59, Tt13 = 170, Tt2.5 = 360, Tt3 = 715, Tt4 = 1400,
Tt5 = 890 deg F; bypass flow 120.202 kg/s, core flow 88.451 kg/s, net
thrust 80.07 kN; Qr = 4.326e7 J/kg; cp = 1005 J/(kg K) through the
compressor and 1089 through the turbine.

As in ``HW_7`` the MATLAB used Rankine temperatures as kelvin; here they
are converted. No fuel flow was given, so f comes from the burner energy
balance with eta_b = 1. The nozzle exit states assume convergent nozzles
with no loss after stations 13 and 5. The data do not close the shaft
power balance: the measured turbine temperature drop supplies only about
70 % of the fan plus compressor work (and gives a turbine efficiency of
0.69), so either Tt5/pt5 were read after a further turbine stage or a
reading is off. The nozzle thrust estimate does agree with the measured
thrust to about 5 %.
"""

from unicodes.propulsion_cycles import Gas, compressor_efficiencies, nozzle, turbine_efficiencies

PSI = 6894.76


def kelvin(deg_f):
    return (deg_f + 459.67) / 1.8


p0 = 101325.0
pt = {k: v * PSI for k, v in {"2": 14.7, "13": 26, "25": 63, "3": 200, "4": 190, "5": 28}.items()}
Tt = {k: kelvin(v) for k, v in {"2": 59, "13": 170, "25": 360, "3": 715, "4": 1400, "5": 890}.items()}
m_fan, m_core, Fn, Qr = 120.202, 88.451, 8.007e4, 4.326e7
c = Gas(1.4, 1005.0)
t = Gas(1089.0 / (1089.0 - 287.0), 1089.0)  # gamma from cp and R = 287

alpha = m_fan / m_core
print(f"bypass ratio alpha = {alpha:.4f}")
for name, a, b in (("fan", "2", "13"), ("LPC (core)", "2", "25"), ("HPC", "25", "3"), ("overall compressor", "2", "3")):
    pi, tau = pt[b] / pt[a], Tt[b] / Tt[a]
    eta, e = compressor_efficiencies(pi, tau, c)
    print(f"{name:20s} pi = {pi:.4f}, tau = {tau:.4f}, eta = {eta:.4f}, e = {e:.4f}")
pi_t, tau_t = pt["5"] / pt["4"], Tt["5"] / Tt["4"]
eta_t, e_t = turbine_efficiencies(pi_t, tau_t, t)
print(f"{'turbine':20s} pi = {pi_t:.4f}, tau = {tau_t:.4f}, eta = {eta_t:.4f}, e = {e_t:.4f}")
f = (t.cp * Tt["4"] - c.cp * Tt["3"]) / (Qr - t.cp * Tt["4"])
print(f"burner: pi_b = {pt['4'] / pt['3']:.4f}, f = {f:.5f} (eta_b = 1), fuel flow {f * m_core:.3f} kg/s")
P_fan = m_fan * c.cp * (Tt["13"] - Tt["2"])
P_comp = m_core * c.cp * (Tt["3"] - Tt["2"])
P_turb = m_core * (1 + f) * t.cp * (Tt["4"] - Tt["5"])
print(f"fan power {P_fan / 1e6:.2f} MW + compressor power {P_comp / 1e6:.2f} MW, turbine power "
      f"{P_turb / 1e6:.2f} MW -> eta_m = {(P_fan + P_comp) / P_turb:.4f}")
e19 = nozzle(pt["13"], Tt["13"], p0, c, convergent=True)
e9 = nozzle(pt["5"], Tt["5"], p0, t, convergent=True)
F_fan = m_fan * (e19.V + e19.pressure_thrust_per_mass(p0))
F_core = m_core * (1 + f) * (e9.V + e9.pressure_thrust_per_mass(p0))
for name, e in (("fan nozzle", e19), ("core nozzle", e9)):
    print(f"{name}: {'choked' if e.choked else 'unchoked'}, M = {e.M:.3f}, V = {e.V:.1f} m/s, p = {e.p / 1e3:.1f} kPa")
print(f"thrust: fan {F_fan / 1e3:.2f} kN + core {F_core / 1e3:.2f} kN = {(F_fan + F_core) / 1e3:.2f} kN "
      f"(measured {Fn / 1e3:.2f})")
print(f"TSFC = {f * m_core / Fn * 1e6:.3f} mg/(N s)")
