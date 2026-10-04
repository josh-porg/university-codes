"""AE 573 final exam problem 1 (``FinalExam_1``): mixed-flow turbofan without afterburner.

M0 = 2.2, p0 = 10 kPa, T0 = 223.15 K, core air flow 25 kg/s; pi_d = 0.85,
pi_c = 15 (overall), pi_f = 1.9 (e_c = e_f = 0.9), pi_b = 0.95
(eta_b = 0.98, Qr = 42.8 MJ/kg), Tt4 = 1673 K, e_t = 0.92, eta_m = 0.95,
core Mach number at the mixer M6 = 0.5, mixer friction loss 0.98,
pi_n = 1, fully expanded nozzle; gamma = 1.4 and cp = 1004 everywhere (as in
the exam). The file also held the student's answers f = 0.013,
mdot_6 = 25.325 kg/s (= 25 (1 + f), so the 25 kg/s is the core flow) and
V9 = 1097.63 m/s; they are printed below for comparison. With these
givens the burner energy balance gives f = 0.0159, not 0.013.

The bypass ratio is not given: it follows from matching the total
pressures at the mixer, ``pt16 = pt6``. Errors in the MATLAB relations:
``a_9^2 = gamma R Tt_9 / (1 + (gamma-1) M_9 / 2)`` and the ideal mixer
pressure ratio used M instead of M^2; ``alpha = (...) / tau_r*(tau_f-1)``
multiplied by ``(tau_f - 1)`` instead of dividing; the ``M_15`` relation
mixed gammas; and typos (``et_m``, ``h_t6``, ``pi_df``) left variables
that nothing else defined. The mixer here is the ideal constant-area
mixer (mass, energy and impulse conserved) times the friction loss.
"""

from unicodes.propulsion_cycles import AIR, mixed_turbofan

m_core = 25.0
r = mixed_turbofan(2.2, 223.15, 10e3, pi_c=15, pi_f=1.9, Tt4=1673, Qr=42.8e6, M6=0.5, gas_c=AIR, gas_t=AIR,
                   pi_d=0.85, pi_b=0.95, pi_m_friction=0.98, pi_n=1.0, e_c=0.9, e_f=0.9, e_t=0.92, eta_b=0.98,
                   eta_m=0.95)
print(r.summary())
m0 = m_core * (1 + r.alpha)
print(f"core flow {m_core:.3f} kg/s, bypass {m0 - m_core:.3f} kg/s, core flow into the mixer "
      f"{m_core * (1 + r.f):.3f} kg/s (file: 25.325), fuel {r.f * m_core:.4f} kg/s")
print(f"f = {r.f:.5f} (file: 0.013), V9 = {r.exits['9'].V:.2f} m/s (file: 1097.63), "
      f"thrust = {r.thrust(m0) / 1e3:.2f} kN")
