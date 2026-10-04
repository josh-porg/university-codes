"""Lift and moment surfaces over (alpha, camber) and the PAH effective curves
(``PAH_Lift_Surface_Gradients``, ``PAH_Moment_Surface_Gradients``, ``PAH_generator``,
``PAH_generation_ploter``, ``airfoil_file_inspector``, ``ReadAirfoil``).

Needs the ``.airfoil`` files (``Selig 1223.airfoil``, ``KL_sexy_sexy.airfoil``)
and, for the tail blanking curve, ``Archytas_horizontal_Tail_Effectiveness.mat``;
pass their folder as the first argument.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from unicodes.aero.airfoil import read_airfoil_file, relaxed_airfoil
from unicodes.aero.pah import CamberSurface, linear_camber_schedule, pitching_moment

folder = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
cambered = read_airfoil_file(folder / "Selig 1223.airfoil")
neutral = read_airfoil_file(folder / "KL_sexy_sexy.airfoil")
relaxed_airfoil(neutral, 1.5).write_airfoil_file(folder / "PAH_KL_sexy_sexy_R1.5.airfoil")

camber_neutral, camber_max = 0.03767, 0.81
surf = CamberSurface.from_airfoils(neutral, cambered, camber_neutral, camber_max, counts=(40, 20, 8))
surf.values["c_l"] /= 1.05  # 3-D correction used in the MATLAB (Polhamus vs 2 pi at A = 30)

try:
    from unicodes.io import load_mat

    d = load_mat(folder / "Archytas_horizontal_Tail_Effectiveness.mat")
    a_eta, eta = np.deg2rad(np.ravel(d["alphas"]))[:-1], np.ravel(d["eta_h"])[:-1]
    eta_h = lambda a: np.interp(a, a_eta, eta) + 0.2  # noqa: E731  (bias correction from the MATLAB)
except FileNotFoundError:
    eta_h = lambda a: 1.0  # noqa: E731

geom = dict(i_h=np.deg2rad(-2.5), C_m0_wf=0.140, C_m0_h=0.0, C_L_alpha_h=0.484, x_cg=0.1217, x_ac_wf=0.185,
            x_ac_h=8.19, S_h=1.75, S=13.27, eta_h=eta_h, downwash_grad=0.111)
cm_n = pitching_moment(surf.alpha[:, 0], neutral, **geom)
cm_c = pitching_moment(surf.alpha[:, -1], cambered, **geom)
w = np.linspace(0, 1, surf.alpha.shape[1])
surf.values["c_m"] = cm_n[:, None] + (cm_c - cm_n)[:, None] * w

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(1, 2, 1, projection="3d")
ax.plot_surface(np.rad2deg(surf.alpha), surf.camber, surf.values["c_l"], cmap="cool", alpha=0.6)
ax.set(xlabel=r"$\alpha$ (deg)", ylabel="camber", zlabel="$c_l$")
alpha = np.linspace(surf.alpha.min(), surf.alpha.max(), 200)
ax2 = fig.add_subplot(1, 2, 2)
for slope in np.linspace(0, -0.15 * 57, 3):
    sched = linear_camber_schedule(0.121, (camber_neutral + camber_max) / 2, slope)
    ax2.plot(np.rad2deg(alpha), surf.slice("c_l", alpha, sched), label=f"dcamber/dalpha = {slope:.2f}")
ax2.set(xlabel=r"$\alpha$ (deg)", ylabel="$c_l$")
ax2.legend()
dCL = surf.slope_along_schedule("c_l", -0.15 * 57)
dCM = surf.slope_along_schedule("c_m", -0.15 * 57)
print("effective C_L_alpha range along the schedule:", np.nanmin(dCL), np.nanmax(dCL))
print("effective C_m_alpha range along the schedule:", np.nanmin(dCM), np.nanmax(dCM))
plt.show()
