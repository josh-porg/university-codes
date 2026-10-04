"""AE 573 homework 8 part 3 (``H_8_part_3``): turbofan fan and compressor, fan nozzle.

M0 = 0.8, p0 = 30 kPa, T0 = 228 K; bypass flow 100 kg/s, core flow 20 kg/s;
pi_d = 0.95, fan pi_f = 1.5, overall compressor pi_c = 32 (e_f = e_c = 0.9),
fan nozzle pi_fn = 0.98 (convergent); gamma = 1.4, R = 287.

The MATLAB mixed ``Tt_1``/``Tt_2`` in the power relations
(``P_fan = mdot_fan cp (Tt_1 - Tt_13)`` as a positive quantity, while the
compressor used ``Tt_2 - Tt_3``), set ``pi_f`` equal to both ``pt_13/pt_1``
and ``pt_2/pt_1`` (putting the whole core behind the fan exit pressure)
and had ``M_19`` use ``gamma_t``. Here the core compressor ratio is the
overall one from the engine face, and powers are positive work inputs.
"""

from unicodes.propulsion_cycles import Gas, compressor, freestream, nozzle

gas = Gas(1.4, 1.4 * 287 / 0.4)
M0, p0, T0 = 0.8, 30e3, 228.0
m_fan, m_core = 100.0, 20.0
pi_d, pi_f, pi_c, e_f, e_c, pi_fn = 0.95, 1.5, 32.0, 0.9, 0.9, 0.98

V0, st0, tau_r, pi_r = freestream(M0, T0, p0, gas)
pt2, Tt2 = st0.pt * pi_d, st0.Tt
tau_f, eta_f = compressor(pi_f, e_f, gas)
tau_c, eta_c = compressor(pi_c, e_c, gas)
Tt13, pt13 = Tt2 * tau_f, pt2 * pi_f
Tt3, pt3 = Tt2 * tau_c, pt2 * pi_c
P_fan = m_fan * gas.cp * (Tt13 - Tt2)
P_comp = m_core * gas.cp * (Tt3 - Tt2)
ex = nozzle(pt13 * pi_fn, Tt13, p0, gas, convergent=True)
F_fan = m_fan * (ex.V - V0 + ex.pressure_thrust_per_mass(p0))

print(f"V0 = {V0:.2f} m/s, Tt0 = Tt2 = {Tt2:.2f} K, pt0 = {st0.pt:.0f} Pa, pt2 = {pt2:.0f} Pa, alpha = {m_fan / m_core:.1f}")
print(f"fan: tau_f = {tau_f:.4f}, eta_f = {eta_f:.4f}, Tt13 = {Tt13:.2f} K, pt13 = {pt13:.0f} Pa, power {P_fan / 1e6:.3f} MW")
print(f"compressor: tau_c = {tau_c:.4f}, eta_c = {eta_c:.4f}, Tt3 = {Tt3:.2f} K, pt3 = {pt3:.0f} Pa, "
      f"power {P_comp / 1e6:.3f} MW")
print(f"total shaft power required {(P_fan + P_comp) / 1e6:.3f} MW")
print(f"fan nozzle: {'choked' if ex.choked else 'unchoked'}, M19 = {ex.M:.4f}, T19 = {ex.T:.2f} K, "
      f"p19 = {ex.p:.0f} Pa, V19 = {ex.V:.2f} m/s")
print(f"fan stream net thrust = {F_fan / 1e3:.2f} kN ({F_fan / m_fan:.2f} N s/kg of bypass air)")
