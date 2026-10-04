"""Jackson Torok, rocket propulsion (``CompiledCode771``, ``JacksonRocketPropullsionCompiledCode771``).

The MATLAB collected the ideal-rocket relations as a symbolic system
(effective exhaust velocity ``c = v2 + (p2 - p3) A2 / mdot``, ``F = c mdot``,
``Isp = c / g0``, mass ratios, ``mdot = (m0 - mf) / t``, ...) and handed
whatever knowns were set to ``solve``. The ``unsorted`` copy set
m0 = 200 kg, mf = 130 kg, payload 110 kg, burn time 3 s, Isp = 240 s; the
other copy set none and then plotted an undefined altitude vector. Both
parts are worked here explicitly:

1. the small motor from the ``unsorted`` knowns (fully expanded, p2 = p3);
2. the commented-out "problem 3": Minuteman first-stage thrust and
   specific impulse versus altitude from the sea-level values (194 600 lbf,
   Isp 254 s, exit pressure 8.66 psia, exit area 1642 in^2), with
   ``F(h) = mdot v2 + (p2 - p_a(h)) A2``.
"""

import matplotlib.pyplot as plt
import numpy as np

# 1. small motor
m0, mf, mpl, t_burn, Isp, g0 = 200.0, 130.0, 110.0, 3.0, 240.0, 9.81
mdot = (m0 - mf) / t_burn
c = Isp * g0
F = c * mdot
print("Part 1:")
print(f"  mdot = {mdot:.4f} kg/s, c = {c:.2f} m/s, F = {F:.1f} N, total impulse = {F * t_burn:.1f} N s")
print(f"  mass ratios: overall mf/m0 = {mf / m0:.4f}, engine (mf - mpl)/(m0 - mpl) = {(mf - mpl) / (m0 - mpl):.4f}, "
      f"propellant fraction zeta = {(m0 - mf) / (m0 - mpl):.4f}")
print(f"  burnout acceleration F/mf = {F / mf:.2f} m/s^2 ({F / mf / g0:.2f} g), ideal delta-v = {c * np.log(m0 / mf):.1f} m/s")

# 2. Minuteman thrust versus altitude (English units, as in the homework)
alt = 3.28084 * np.array([0, 1000, 3000, 5000, 10000, 25000, 50000, 75000, 100000, 130000, 160000, 200000, 300000,
                          400000, 600000, 1000000])
pr = np.array([1, 0.887, 0.66919, 0.53313, 0.26151, 0.025158, 0.00078735, 2.0408e-5, 3.15983e-7, 1.2341e-8,
               2.9997e-9, 8.3628e-10, 8.6557e-11, 1.4328e-11, 8.1056e-13, 7.4155e-14])
p_sl = 101325 * 0.020885  # psf
p2, F_sl, A2, Isp_sl, g0e = 8.66 * 144, 194600.0, 1642 / 144, 254.0, 32.174
mdot = F_sl / (Isp_sl * g0e)  # slug/s
v2 = (F_sl - (p2 - p_sl) * A2) / mdot
F_alt = mdot * v2 + (p2 - p_sl * pr) * A2
Isp_alt = F_alt / (mdot * g0e)
print("Part 2:")
print(f"  mdot = {mdot:.2f} slug/s, exit velocity {v2:.1f} ft/s; vacuum thrust {F_alt[-1]:.0f} lbf, Isp {Isp_alt[-1]:.1f} s")
fig, ax1 = plt.subplots()
ax1.plot(alt, F_alt / 1000, color=(0.75, 0.25, 0))
ax1.set(ylabel="Thrust (lb x10^3)", xlabel="Altitude (ft)", ylim=(190, 220), title="Minuteman Rocket Performance")
ax2 = ax1.twinx()
ax2.plot(alt, Isp_alt, color=(0, 0.45, 0.45))
ax2.set(ylabel="Specific Impulse (sec)", ylim=(250, 290))
ax1.set_xscale("log")
plt.show()
