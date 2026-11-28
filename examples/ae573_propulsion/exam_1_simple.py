"""AE 573 exam 1 (``EXAM_1_simple``): separate-exhaust turbofan.

M0 = 0.86, p0 = 30 kPa, T0 = 248 K; bypass 240 kg/s, core 30 kg/s
(alpha = 8); pi_d = 0.99, pi_f = 1.4, pi_c = 45 (overall, e_f = e_c = 0.9),
pi_b = 0.97 (eta_b = 0.99, Qr = 42 MJ/kg), Tt4 = 1625 K, e_t = 0.8,
eta_m = 0.99, pi_n = pi_fn = 0.98, convergent nozzles (the fan nozzle was
given as choked, M19 = 1).

Errors in the MATLAB file: ``gamma`` was assigned 1.4 and then 1.33, so the
whole engine (inlet and fan included) ran with gamma = 1.33; ``e_t = 1/.8``
gave a turbine polytropic efficiency of 1.25, which, with the turbine
relation written as ``tau_t = pi_t^((gamma-1)/gamma * e_t)``, had the
turbine producing temperature rise (Tt5 = 2416 K > Tt4); and the
fuel heating value appeared twice in different units (``Qr_4 = 42000``). Here
the cold sections use gamma = 1.4, cp = 1004 and the hot sections gamma =
1.33, cp = 1156, with e_t = 0.8.
"""

from unicodes.propulsion_cycles import turbofan

m_fan, m_core = 240.0, 30.0
r = turbofan(0.86, 248.0, 30e3, pi_c=45, pi_f=1.4, alpha=m_fan / m_core, Tt4=1625, Qr=42e6, pi_d=0.99, pi_b=0.97,
             pi_n=0.98, pi_fn=0.98, e_c=0.9, e_f=0.9, e_t=0.8, eta_b=0.99, eta_m=0.99, convergent=True)
print(r.summary())
F = r.thrust(m_fan + m_core)
print(f"thrust = {F / 1e3:.2f} kN, fuel flow = {r.f * m_core:.4f} kg/s")
