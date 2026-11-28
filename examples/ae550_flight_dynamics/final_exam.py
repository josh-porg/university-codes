"""AE 550 final exam: canard fighter kinematics, lift build-up, trim and lateral derivatives (``final_exam``)."""

import numpy as np

from unicodes.aero import stability, wing
from unicodes.flight_dynamics import kinematics as k

d = np.rad2deg
V_I = np.array([850.0, -30, -30])  # ft/s, NED
att = np.deg2rad([-5, 4, 0])
rho, a, W, CLmax, e, f = 0.00175549, 1077.39, 42900.0, 1.1, 0.92, 10.0
S, A, b, taper, sweep = 492.0, 2.7, 36.75, 0.27, np.deg2rad(37)
M = np.linalg.norm(V_I) / a

gamma, chi = k.flight_path_angles(V_I)
print(f"Q1 gamma {d(gamma):.3f} deg, heading {d(chi):.3f} deg")
H = k.dcm_inertial_to_body(att)
V_B = H @ V_I
alpha, beta = k.air_flow_angles(V_B)
print(f"Q2 V_B = {V_B}; Q3 alpha {d(alpha):.3f}, beta {d(beta):.3f} deg")
for name, wind in [("SE", [80 * np.cos(np.pi / 4), -80 * np.cos(np.pi / 4), 0]), ("NW", [-80 * np.cos(np.pi / 4), 80 * np.cos(np.pi / 4), 0])]:
    Va = V_B - H @ np.array(wind)
    aw, bw = k.air_flow_angles(Va)
    print(f"Q5/8 {name} wind: d alpha {d(aw - alpha):.3f}, d beta {d(bw - beta):.3f} deg, d qbar {0.5 * rho * (Va @ Va - V_I @ V_I):.1f} psf")
print("Q9 max tailwind before stall:", V_B[0] - np.sqrt(2 * W / (rho * CLmax * S)), "ft/s")

C_L_a_wf = wing.polhamus_lift_slope(2 * np.pi, A, sweep, taper, M)
C_L_a_c = wing.polhamus_lift_slope(2 * np.pi, 1.0, np.deg2rad(55), 0.19, M)
C_L_a_v = wing.polhamus_lift_slope(2 * np.pi, 3.2, sweep, 0.30, M)
ratio = C_L_a_wf / wing.polhamus_lift_slope(2 * np.pi, A, sweep, taper, 0.0)
upwash = wing.downwash_gradient(0.095, 0.9176, A, sweep, taper, ratio)
print(f"Q10 C_L_alpha wf {C_L_a_wf:.4f}, canard {C_L_a_c:.4f}, fin {C_L_a_v:.4f}; Q11 de/da = {upwash:.4f}")

canard = stability.Surface(S=93.7, C_L_alpha=C_L_a_c, x_ac=-0.8571, eta=0.98, downwash_grad=upwash, canard=True)
lon = stability.longitudinal(S, 0.08, C_L_a_wf, 0.25, 0.025, x_cg=0.0, surface=canard)
x_cg = lon.x_ac - 0.05
lon = stability.longitudinal(S, 0.08, C_L_a_wf, 0.25, 0.025, x_cg=x_cg, surface=canard)
print(f"Q12 C_L_alpha {lon.C_L_alpha:.4f}, C_L_ic {lon.C_L_i:.4f}; Q16 x_ac {lon.x_ac:.4f}; Q17 x_cg {x_cg:.4f}")
print(f"Q18 C_m0 {lon.C_m0:.4f}, C_m_alpha {lon.C_m_alpha:.4f}, C_m_ic {lon.C_m_i:.4f}; Q19 i_c trim {d(lon.trim_incidence(alpha)):.3f} deg")
bad = stability.longitudinal(S, 0.08, C_L_a_wf, 0.25, 0.025, x_cg=0.55, surface=canard)
print(f"Q20 cg at 0.55: static margin {bad.static_margin:.3f}, trim {d(bad.trim_incidence(alpha)):.3f} deg")

C_D0 = f / S
for theta in (0, 15, 30):
    Vb = k.dcm_inertial_to_body(np.deg2rad([-5, theta, 0])) @ V_I
    al, _ = k.air_flow_angles(Vb)
    C_L = lon.C_L0 + lon.C_L_alpha * al
    T = stability.thrust_required(0.5 * rho * Vb @ Vb, S, C_L, C_D0 + C_L**2 / (np.pi * A * e), al, W, np.deg2rad(theta))
    print(f"Q14 theta {theta:2d}: alpha {d(al):6.2f} deg, T/W = {T / W:.3f}")

lat = stability.lateral(S, b, (-0.6, 0.05, 0.04),
                        dict(S=88.3, C_L_alpha=C_L_a_v, eta=0.93, sidewash_grad=0.18, x=14.256, z=6.48, tau_r=0.5),
                        horizontal=(93.7, 21.315, 0.98, -0.12, -0.01, 0.01))
print("Q21", lat)
print(f"Q23 rudder {d(lat.rudder_for_sideslip(beta)):.3f} deg, aileron {d(lat.aileron_for_sideslip(beta, 0.05)):.3f} deg")
