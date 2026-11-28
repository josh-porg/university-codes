"""AE 573 turbojet homework (``HomeWork``): performance from a measured fuel flow and exhaust velocity.

mdot_0 = 100 kg/s at M0 = 2, p0 = 20 kPa, T0 = 228 K; fuel 2.8 kg/s
(Qr = 42 MJ/kg, eta_b = 0.995), fully expanded exhaust at V9 = 1200 m/s;
pi_d = 0.9, pi_c = 12 (e_c = 0.9), pi_b = 0.95; cold gas gamma 1.4,
cp 1004, hot gas gamma 1.33, cp 1156.

The MATLAB fed ~110 overlapping relations to the relation solver. Several
were wrong (``pi_d = pt_0 / pt_2`` is inverted; ``tau_lambda = cp_4*Tt_4 /
cp_0*T_0`` multiplies by T_0; the afterburner balance ``f_AB = (1+f)(ht_7 -
ht_5) / (Qr_AB eta_AB - ht_7)`` with ``Qr_AB = 0``; ``eta_p = 2 / (1 + V_9 /
V_0)`` only holds for f = 0) and the result depended on which relation the
solver reached first. Here the quantities follow in order from the givens.
The burner exit temperature comes from the fuel flow through the burner
energy balance and the turbine exit temperature from the shaft power
balance (eta_m = 1).
"""

from unicodes.propulsion_cycles import AIR, COMBUSTION_GAS, compressor, freestream

mdot0, M0, p0, T0 = 100.0, 2.0, 20e3, 228.0
mdot_f, Qr, eta_b, V9 = 2.8, 42e6, 0.995, 1200.0
pi_d, pi_c, e_c, pi_b = 0.9, 12.0, 0.9, 0.95
c, t = AIR, COMBUSTION_GAS
g0 = 9.81

V0, st0, tau_r, pi_r = freestream(M0, T0, p0, c)
f = mdot_f / mdot0
Fn = (mdot0 + mdot_f) * V9 - mdot0 * V0  # p9 = p0
Dram = mdot0 * V0
TSFC = mdot_f / Fn
dKE = ((mdot0 + mdot_f) * V9**2 - mdot0 * V0**2) / 2
eta_th = dKE / (mdot_f * Qr)
eta_p = Fn * V0 / dKE
pt2, Tt2 = st0.pt * pi_d, st0.Tt
tau_c, eta_c = compressor(pi_c, e_c, c)
Tt3, pt3 = Tt2 * tau_c, pt2 * pi_c
Tt4 = (eta_b * f * Qr + c.cp * Tt3) / ((1 + f) * t.cp)  # (1+f) cp_t Tt4 - cp_c Tt3 = eta_b f Qr
Tt5 = Tt4 - c.cp * (Tt3 - Tt2) / ((1 + f) * t.cp)
P_compressor = mdot0 * c.cp * (Tt3 - Tt2)

print(f"V0 = {V0:.2f} m/s, tau_r = {tau_r:.4f}, pi_r = {pi_r:.4f}, Tt0 = {st0.Tt:.2f} K, pt0 = {st0.pt:.0f} Pa")
print(f"ram drag = {Dram / 1e3:.2f} kN, gross thrust = {(mdot0 + mdot_f) * V9 / 1e3:.2f} kN, Fn = {Fn / 1e3:.2f} kN")
print(f"F/mdot0 = {Fn / mdot0:.2f} N s/kg, TSFC = {TSFC * 1e6:.3f} mg/(N s), Isp = {Fn / (mdot_f * g0):.1f} s")
print(f"eta_th = {eta_th:.4f}, eta_p = {eta_p:.4f}, eta_o = {eta_th * eta_p:.4f}")
print(f"pt2 = {pt2:.0f} Pa, Tt2 = {Tt2:.2f} K; compressor tau_c = {tau_c:.4f}, eta_c = {eta_c:.4f}")
print(f"Tt3 = {Tt3:.2f} K, pt3 = {pt3:.0f} Pa, compressor power = {P_compressor / 1e6:.2f} MW")
print(f"Tt4 = {Tt4:.2f} K (from the fuel flow), pt4 = {pt3 * pi_b:.0f} Pa, Tt5 = {Tt5:.2f} K")
