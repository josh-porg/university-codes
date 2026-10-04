"""AE 573 practice quiz (``Pratice_quiz_turbofan``): separate-exhaust turbofan at cruise.

M0 = 0.8 at 11.2776 km (37 000 ft, standard atmosphere); alpha = 6;
pi_d = 0.995, pi_f = 1.6, pi_c = 40 (e_f = e_c = 0.9), pi_b = 0.95
(eta_b = 0.98, Qr = 42.8 MJ/kg), tau_lambda = cp_t Tt4 / (cp_c T0) = 7,
e_t = 0.9, eta_m = 0.975, pi_n = pi_fn = 0.98, convergent nozzles choked
(M9 = M19 = 1); gamma 1.4 / 1.33, cp 1004 / 1156 J/(kg K).

The MATLAB gave cp in kJ/(kg K) and Qr in kJ/kg but R = 287 J/(kg K),
so every relation joining cp and R (speeds of sound, ``gamma = cp / cv``)
mixed units. That is where the negative ``u_19 = -322 422`` came from. Its
``tau_lambda = cp_4*Tt_4 / cp_0*T_0`` also multiplied by T0 instead of
dividing. Here everything is in SI.
"""

from unicodes.atmosphere import isa
from unicodes.propulsion_cycles import AIR, COMBUSTION_GAS, turbofan

atm = isa(11.2776e3)
T0, p0 = atm.T, atm.p
Tt4 = 7 * AIR.cp * T0 / COMBUSTION_GAS.cp
print(f"T0 = {T0:.2f} K, p0 = {p0:.0f} Pa, Tt4 = {Tt4:.1f} K (tau_lambda = 7)")
r = turbofan(0.8, T0, p0, pi_c=40, pi_f=1.6, alpha=6, Tt4=Tt4, Qr=42.8e6, pi_d=0.995, pi_b=0.95, pi_n=0.98,
             pi_fn=0.98, e_c=0.9, e_f=0.9, e_t=0.9, eta_b=0.98, eta_m=0.975, convergent=True)
print(r.summary())
