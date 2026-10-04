"""AE 550 final exam part 2: conventional-tail business jet stability and trim (``final_exam_part_2``)."""

import numpy as np

from unicodes.aero import stability

d = np.rad2deg
a, M, rho, W = 1097.09, 0.9, 0.00204, 11700.0
alpha, beta = np.deg2rad(0.5), np.deg2rad(2)
S, b, c_bar, A, e = 170.0, 25.3, 7.4, 25.3**2 / 170, 0.9
ct_h, cr_h, b_h = 1.8512, 6.1707, 13.8841
S_h = b_h / 2 * cr_h * (1 + ct_h / cr_h)
ct_v, cr_v, b_v = 2.4683, 9.5646, 6.7878
S_v = b_v / 2 * cr_v * (1 + ct_v / cr_v)

tail = stability.Surface(S=S_h, C_L_alpha=5.82, x_ac=12.3415 / c_bar, eta=0.9, downwash_grad=0.3907, epsilon_0=np.deg2rad(0.4))
lon = stability.longitudinal(S, 0.12, 3.1, 0.25, -0.08, x_cg=0.6, surface=tail, ac_shift=-0.08)
print(f"Q1 C_L0 {lon.C_L0:.4f}, C_L_alpha {lon.C_L_alpha:.4f}, C_L_ih {lon.C_L_i:.4f}")
C_L = lon.C_L0 + lon.C_L_alpha * alpha
C_D = 0.0162 + C_L**2 / (np.pi * A * e)
q = 0.5 * rho * (M * a) ** 2
print(f"Q2 level thrust {stability.thrust_required(q, S, C_L, C_D, alpha, W, 0.0):.0f} lbf")
print(f"Q3 C_m0 {lon.C_m0:.4f}, C_m_alpha {lon.C_m_alpha:.4f}, C_m_ih {lon.C_m_i:.4f}")
print(f"Q4 neutral point {lon.x_ac:.4f}, static margin {lon.static_margin:.4f}")
moved = stability.longitudinal(S, 0.12, 3.1, 0.4, -0.08, x_cg=0.6,
                               surface=stability.Surface(S_h, 5.82, 12.0329 / c_bar, 0.9, 0.3907, np.deg2rad(0.4)), ac_shift=-0.08)
print(f"Q5 moved a.c.: static margin {moved.static_margin:.4f}; Q7 trim i_h {d(lon.trim_incidence(alpha)):.3f} deg")
lat = stability.lateral(S, b, (-0.23, -0.016, 0.2),
                        dict(S=S_v, C_L_alpha=1.73, eta=0.96, sidewash_grad=0.3, x=10.7988, z=4.6280, tau_r=0.4),
                        horizontal=(S_h, b_h, 0.9, -0.04, -0.009, 0.03))
print("Q8", lat)
print(f"Q9 rudder {d(lat.rudder_for_sideslip(beta)):.3f} deg; Q10 aileron {d(lat.aileron_for_sideslip(beta, 0.013)):.3f} deg")
