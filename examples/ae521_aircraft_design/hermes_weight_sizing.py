"""Roskam fuel-fraction weight sizing for the Hermes-launch mission (``vehicle_weight_sizing_1``, ``weight_Sizing_Hermes_V_0``).

Mission 1: climb to 65,000 ft, release two Hermes rockets, loiter 8 min,
descend, balked landing, divert 100 nmi, loiter 45 min, land. British
units as in the original.
"""

import numpy as np

from unicodes.aero.sizing import breguet_endurance_fraction, breguet_range_fraction, estimate_weights, mach_to_kts

L_D, c_j = 26.0, 0.34
payload = (30000 + 30000 * 0.25) * 2  # two Hermes rockets
crew = 175 * 4
mission = [
    ("engine start", 0.99), ("taxi", 0.99), ("take-off", 0.995), ("climb", 0.98),
    ("loiter", breguet_endurance_fraction(8 / 60, c_j, L_D)),
    ("descent", 0.995), ("balked landing", 1.0), ("short climb", 0.995),
    ("divert", breguet_range_fraction(100, mach_to_kts(0.75, 1036.8), c_j, L_D)),
    ("loiter 45 min", breguet_endurance_fraction(0.75, c_j, L_D)),
    ("descent 2", 0.995), ("landing", 1.0), ("taxi and shutdown", 0.995),
]
est = estimate_weights(140000, mission, payload=payload, crew=crew, payload_drop_after="climb")
print(est)

# weight_Sizing_Hermes_V_0: Hermes launch weight for a 200 nmi Mach 4 cruise at L/D 4, Isp 250 s
a = 968.1  # ft/s at 65,000 ft
V, I_sp, R = 4 * a, 250.0, 200 * 6076.0
W_p = (30000 + 30000) / 2
k = np.exp(R / (V * 4 * I_sp))
W_L = -W_p * k / (1 - k)
print(f"Hermes launch weight {W_L:.0f} lbf, propellant fraction {W_p / W_L:.3f}, dry weight {W_L - W_p:.0f} lbf")
