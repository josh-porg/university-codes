"""AE 573 homework 7 (``HW_7``): component performance from turbojet test-stand data.

Static sea-level test (V0 = 0, p0 = 101.325 kPa) of a two-spool turbojet
in English units: pt2 = 14.7, pt2.5 = 54, pt3 = 167, pt4 = 158, pt5 = 36,
pt9 = 33 psia; Tt2 = 59, Tt2.5 = 330, Tt3 = 660, Tt4 = 1570,
Tt5 = Tt9 = 1013 deg F; fuel 8520 lb/h (Qr = 4.326e7 J/kg), air 164 lb/s,
thrust 10 200 lbf.

The MATLAB converted the temperatures to Rankine (F + 460) and then used
them as kelvin, so every temperature-based result (the burner balance,
speeds of sound, exit velocity) was off by a factor of 1.8. It also had
``tau_c_lp = Tt_2 / Tt_2_5`` upside down and set ``A_9 = 0``. Here the
temperatures are converted to kelvin; cold sections use gamma = 1.4,
cp = 1004 and hot sections gamma = 1.33, cp = 1156.
"""

from unicodes.propulsion_cycles import (
    AIR,
    COMBUSTION_GAS,
    compressor_efficiencies,
    nozzle,
    turbine_efficiencies,
)
from unicodes.thermo.combustion import enthalpy_molar

PSI, LBM, LBF = 6894.76, 0.453592, 4.44822


def kelvin(deg_f):
    return (deg_f + 459.67) / 1.8


c, t = AIR, COMBUSTION_GAS
p0 = 101325.0
pt = {k: v * PSI for k, v in {"2": 14.7, "25": 54, "3": 167, "4": 158, "5": 36, "9": 33}.items()}
Tt = {k: kelvin(v) for k, v in {"2": 59, "25": 330, "3": 660, "4": 1570, "5": 1013, "9": 1013}.items()}
mdot_f, mdot0, F = 8520 * LBM / 3600, 164 * LBM, 10.2e3 * LBF
Qr = 4.326e7

f = mdot_f / mdot0
print(f"mdot0 = {mdot0:.2f} kg/s, mdot_f = {mdot_f:.4f} kg/s, f = {f:.5f}, F = {F / 1e3:.2f} kN")
for name, a, b in (("LPC", "2", "25"), ("HPC", "25", "3"), ("overall compressor", "2", "3")):
    pi, tau = pt[b] / pt[a], Tt[b] / Tt[a]
    eta, e = compressor_efficiencies(pi, tau, c)
    print(f"{name:20s} pi = {pi:.4f}, tau = {tau:.4f}, eta = {eta:.4f}, e = {e:.4f}")
pi_t, tau_t = pt["5"] / pt["4"], Tt["5"] / Tt["4"]
eta_t, e_t = turbine_efficiencies(pi_t, tau_t, t)
print(f"{'turbine':20s} pi = {pi_t:.4f}, tau = {tau_t:.4f}, eta = {eta_t:.4f}, e = {e_t:.4f}")
pi_b = pt["4"] / pt["3"]
eta_b = ((1 + f) * t.cp * Tt["4"] - c.cp * Tt["3"]) / (f * Qr)
# The two-constant-cp model jumps from cp_c at Tt3 to cp_t at Tt4 and overstates the heat added;
# with the real enthalpy of air (NASA polynomials) between the same temperatures it is consistent.
M_air = 0.79 * 28.013 + 0.21 * 31.999
h_air = lambda T: (0.79 * enthalpy_molar("N2", T) + 0.21 * enthalpy_molar("O2", T)) / M_air * 1e3  # noqa: E731
eta_b_real = (1 + f) * (h_air(Tt["4"]) - h_air(Tt["3"])) / (f * Qr)
print(f"burner: pi_b = {pi_b:.4f}, eta_b = {eta_b:.4f} with cp_c / cp_t (above 1: the constant-cp model "
      f"overstates the heat), {eta_b_real:.4f} with the enthalpy of air")
P_c = mdot0 * c.cp * (Tt["3"] - Tt["2"])
P_t = mdot0 * (1 + f) * t.cp * (Tt["4"] - Tt["5"])
print(f"compressor power {P_c / 1e6:.2f} MW, turbine power {P_t / 1e6:.2f} MW, eta_m = {P_c / P_t:.4f}")
print(f"nozzle pi_n = pt9 / pt5 = {pt['9'] / pt['5']:.4f}")
ex = nozzle(pt["9"], Tt["9"], p0, t, convergent=True)
F_calc = mdot0 * (1 + f) * ex.V + mdot0 * (1 + f) * ex.pressure_thrust_per_mass(p0)
print(f"convergent nozzle: {'choked' if ex.choked else 'unchoked'}, M9 = {ex.M:.3f}, V9 = {ex.V:.1f} m/s, "
      f"p9 = {ex.p / 1e3:.1f} kPa -> F = {F_calc / 1e3:.2f} kN (measured {F / 1e3:.2f})")
V_eff = F / (mdot0 * (1 + f))
dKE = mdot0 * (1 + f) * V_eff**2 / 2
print(f"effective exhaust velocity F / mdot9 = {V_eff:.1f} m/s, TSFC = {mdot_f / F * 1e6:.3f} mg/(N s), "
      f"Isp = {F / (mdot_f * 9.80665):.0f} s, eta_th = {dKE / (mdot_f * Qr):.4f}")
print(f"Carnot limit between Tt2 and Tt4: {1 - Tt['2'] / Tt['4']:.4f}")
