"""AE 508 final exam problems 3, 4 and 5: FE stress-gradient extrapolation and mesh convergence.

The FE results (NASTRAN/Patran output typed into the MATLAB scripts) are
embedded. Stresses near a singular boundary are extrapolated to the edge by
fitting the smooth part of the gradient (MATLAB ``fit`` -> ``np.polyfit``).
"""

import matplotlib.pyplot as plt
import numpy as np

# Problem 3: minimum principal stress gradient for three meshes (ConvergencePlots)
fine = np.array([-3.103425e5, -3.103011e5, -3.102348e5, -3.101149e5, -3.100142e5, -3.102999e5, -3.114438e5,
                 -3.123515e5, -3.073643e5, -2.891331e5, -2.744139e5, -3.322485e5])
medium = np.array([-3.063436e5, -3.058102e5, -3.064954e5, -3.093698e5, -2.897105e5, -3.191785e5])
coarse = np.array([-2.900051e5, -2.861737e5, -2.928853e5, -2.900515e5])
d_fine = np.arange(11, -1, -1) / 4
d_medium = d_fine[1::2]
d_coarse = np.r_[3, d_medium[1::2]]
fits = [np.polyfit(d_fine[3:8], fine[3:8], 1), np.polyfit(d_medium[:4], medium[:4], 1),
        np.polyfit(d_medium[:3], coarse[:3], 1)]  # the MATLAB fitted the coarse values at the medium stations
dof = np.array([9616, 2504, 676])
edge = [np.polyval(c, 0) for c in fits]
fig, ax = plt.subplots()
for d, s, c, name in zip((d_fine, d_medium, d_coarse), (fine, medium, coarse), fits, ("fine", "medium", "coarse")):
    line, = ax.plot(d, s, "x", label=name)
    ax.plot(d_fine, np.polyval(c, d_fine), color=line.get_color(), label=f"{name} curve fit")
ax.set(xlabel="Distance from boundary condition (in)", ylabel="Min principal stress (psi)", title="Min principal gradient")
ax.legend()
print("Problem 3 extrapolated edge stress (fine, medium, coarse):", np.round(edge))

maxp = np.array([3.038211e5, 3.273242e5, 3.245280e5])
disp = np.linalg.norm([[1.236948e-3, 6.831180, 8.887047], [4.979953e-9, 7.119062, 9.116643],
                       [2.506017e-9, 7.197069, 9.176082]], axis=1)
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, y, label in zip(axes, (edge, disp, maxp), ("Min principal stress (psi)", "Max displacement (in)",
                                                     "Max principal stress (psi)")):
    dofs = dof if y is edge else dof[::-1]
    ax.semilogx(dofs, y, "xr-")
    ax.set(xlabel="Degrees of freedom", ylabel=label)
fig.suptitle("Problem 3 convergence")

# Problem 4: stress gradient at the lower side of the hole
s4 = np.array([5.752410e4, 4.910493e4, 4.297333e4, 3.831448e4, 3.472026e4, 3.220662e4, 3.051528e4, 2.952276e4,
               2.897803e4, 2.880332e4, 2.896937e4])
d4 = 0.04 * np.arange(11)
p4, p2 = np.polyfit(d4[1:], s4[1:], 4), np.polyfit(d4[2:], s4[2:], 2)
fig, ax = plt.subplots()
ax.plot(d4, s4, "xr", label="FEM data")
dd = np.linspace(0, d4[-1], 200)
ax.plot(dd, np.polyval(p4, dd), "g", label="4th order polynomial")
ax.plot(dd, np.polyval(p2, dd), "b", label="2nd order polynomial")
ax.set(title="Stress gradient of lower side of hole", xlabel="Distance (in)", ylabel="Max principal stress (psi)")
ax.legend()
print(f"Problem 4 edge stress: quartic {np.polyval(p4, 0):.0f} psi, quadratic {np.polyval(p2, 0):.0f} psi")

# Problem 5: two gradient extrapolations and the principal shear stresses
d5 = 0.1 * np.arange(9)
a = np.array([-3.737590e3, -3.659510e3, -3.238332e3, -3.034609e3, -2.901769e3, -2.804090e3, -2.726157e3, -2.661634e3,
              -2.606965e3])
b = np.array([-3.255075e3, -3.561667e3, -3.010950e3, -2.758737e3, -2.509571e3, -2.293920e3, -2.119923e3, -1.979570e3,
              -1.867243e3])
print(f"Problem 5 edge stress: boundary condition {np.polyval(np.polyfit(d5[3:], a[3:], 1), 0):.1f} psi, "
      f"free edge {np.polyval(np.polyfit(d5[1:], b[1:], 2), 0):.1f} psi")
prince = np.array([[5.640320e3, -3.189271e3], [1.950545e3, -1.799015e3], [2.055740e3, -1.125872e3],
                   [1.858906e3, -9.842408e2], [1.726500e3, -9.543512e2], [1.578808e3, -9.777998e2],
                   [1.443706e3, -1.018443e3], [1.318383e3, -1.066849e3], [1.204071e3, -1.117685e3]])
tau12 = (prince[:, 0] - prince[:, 1]) / 2
print("Problem 5 max shear (psi):", np.round(np.max(np.abs([tau12, prince[:, 1] / 2, prince[:, 0] / 2]), axis=0)))
plt.show()
