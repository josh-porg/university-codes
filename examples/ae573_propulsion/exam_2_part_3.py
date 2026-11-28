"""AE 573 exam 2 part 3 (``Exam_2_part_3``): turboprop.

M0 = 0.76, p0 = 35 kPa, T0 = 242 K; pi_d = 0.98, pi_c = 25 (e_c = 0.9),
pi_b = 0.95 (eta_b = 0.99, Qr = 42 MJ/kg), Tt4 = 1600 K; high-pressure
turbine e_hpt = 0.85, eta_m_hpt = 0.99; of the expansion work available
after the HPT, the fraction alpha = 0.87 goes to the free (low-pressure)
turbine (eta_lpt = 0.9, eta_m_lpt = 0.985), which drives the propeller
through a gearbox (eta_gb = 0.996, eta_prop = 0.75); the rest expands in the
nozzle (eta_n = 0.96). Gamma 1.4 / 1.33, cp 1004 / 1156.

Errors in the MATLAB relations: ``(p_0/pt_45)^(gamma-1)/gamma`` lacked
the brackets around the exponent; ``v_9 = 2 (1-alpha) eta_n ...`` was
V9^2, not V9; ``eta_gb = P_prop/P_lpt*eta_gb`` and
``P_lpt = P_prop*eta_prop*eta_gb`` contradicted ``P_prop = eta_gb eta_m_lpt P_lpt``;
``TSFC = mdot_f/P_prop+P_core`` dropped the brackets; ``Qr = 42000`` was in
kJ/kg while cp was in J/(kg K); and the fan relations (``pi_f = 1.5``) do
not belong to a turboprop.
"""

from unicodes.propulsion_cycles import turboprop

r = turboprop(0.76, 242.0, 35e3, pi_c=25, Tt4=1600, Qr=42e6, work_split=0.87, pi_d=0.98, pi_b=0.95, e_c=0.9,
              e_hpt=0.85, eta_b=0.99, eta_m_hpt=0.99, eta_lpt=0.9, eta_m_lpt=0.985, eta_gb=0.996, eta_prop=0.75,
              eta_n=0.96)
print(r.summary())
x = r.extra
print(f"propeller shaft power {x['P_prop/mdot0'] / 1e3:.2f} kW per kg/s of air, propeller thrust "
      f"{x['F_prop/mdot0']:.2f} N s/kg, jet thrust {x['F_core/mdot0']:.2f} N s/kg "
      f"({100 * x['F_core/mdot0'] / r.specific_thrust:.1f} % of the total)")
print(f"power specific fuel consumption = {x['power_specific_fuel_consumption'] * 3.6e9:.1f} g/(kW h)")
