"""Jackson Torok, AE 211 homework 10 (``Torok_Jackson_HW10_MATLAB``): interpolation, curve fitting,
numerical differentiation and integration, and an ODE (13.4-13.20).
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline, PchipInterpolator, interp1d

data = np.array([[1, 2494, 4157], [2, 1247, 2078], [3, 831, 1386], [4, 623, 1039], [5, 499, 831], [6, 416, 693]])
vol, pres = data[:, 0], data[:, 1:]
# 13.4 (MATLAB 'cubic' is shape-preserving pchip; a cubic spline is shown too)

print("13.4: linear at V = 5.2:", interp1d(vol, pres, axis=0)(5.2), " pchip:", PchipInterpolator(vol, pres)(5.2),
      " spline:", CubicSpline(vol, pres)(5.2))
# 13.7
plt.figure(num="13.7")
plt.plot(vol, pres[:, 0], "o")
x = np.arange(1, 6.01, 0.2)
for deg in (1, 2, 3, 4):
    plt.plot(x, np.polyval(np.polyfit(vol, pres[:, 0], deg), x), label=f"degree {deg}")
plt.legend()
# 13.16
time = np.arange(25)
alt = np.array([0, 107.37, 210, 307.63, 400, 484.6, 550, 583.97, 580, 549.53, 570, 699.18, 850, 927.51, 950, 954.51,
                940, 910.68, 930, 1041.52, 1150, 1158.24, 1100, 1041.76, 1050])
fig, ax = plt.subplots(3, 1, num="13.16")
ax[0].plot(time, alt)
ax[0].set_title("Altitude")
ax[1].plot(time[1:], np.diff(alt))
ax[1].set_title("Velocity")
ax[2].plot(time[2:], np.diff(alt, 2))
ax[2].set_title("Acceleration")
print("13.16: acceleration changes sign near t = 8, 13 and 21 s")
# 13.18
a_, b_, c_, d_ = 25.48, 1.522e-2, -0.7155e-5, 1.312e-9
dh = quad(lambda T: a_ + b_ * T + c_ * T**2 + d_ * T**3, 300, 1000)[0]
print(f"13.18: delta h = {dh:.2f} kJ/kmol")
# 13.20: y' = 1 - sin t, y(0) = 1  ->  y = t + cos t
print(f"13.20: y(4) - y(0) = {4 + np.cos(4) - 1:.6f}")
plt.show()
