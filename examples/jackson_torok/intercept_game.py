"""Jackson Torok, AE 211 app (``Untitled``, an App Designer GUI): intercept an incoming ballistic shot.

The enemy round leaves x = 0 at 45 deg and 1735 m/s towards the city at
x = 307 km. Until t = 60 s the player sets three launch angles; at 60 s
three interceptors leave the city at twice the enemy speed. A hit is a
pass within 5 m; the game is lost when the enemy round lands. The GUI
animation is replaced by a fixed-step simulation (1 ms steps within 1 km
of an interceptor, as the app did, 10 ms otherwise) and a plot of the
trajectories; ``--solve`` finds the launch angle that hits (the app's
comment gave 28.643 deg)::

    python intercept_game.py --angles 20 28.643 35
    python intercept_game.py --solve
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar

G, V_ENEMY, LAUNCH_ANGLE, CITY, LAUNCH_TIME = 9.8, 1735.0, 45.0, 307e3, 60.0
V_SHOT = 2 * V_ENEMY


def enemy(t):
    return t * V_ENEMY * np.cos(np.radians(LAUNCH_ANGLE)), t * V_ENEMY * np.sin(np.radians(LAUNCH_ANGLE)) - G * t**2 / 2


def shot(t, angle):
    s = t - LAUNCH_TIME
    return CITY - s * V_SHOT * np.cos(np.radians(angle)), s * V_SHOT * np.sin(np.radians(angle)) - G * s**2 / 2


def play(angles, hit_radius=5.0):
    """Returns ``(hit, time, closest_approach_per_shot)``."""
    t, closest = 0.0, np.full(len(angles), np.inf)
    while True:
        xm, ym = enemy(t)
        if ym < -0.01:
            return False, t, closest
        if t >= LAUNCH_TIME:
            d = np.array([np.hypot(xm - sx, ym - sy) for sx, sy in (shot(t, a) for a in angles)])
            closest = np.minimum(closest, d)
            if np.any(d <= hit_radius):
                return True, t, closest
            t += 1e-3 if np.any(d <= 1000) else 1e-2
        else:
            t += 1e-2


def miss_distance(angle):
    """Closest approach of one interceptor (continuous time)."""
    res = minimize_scalar(lambda t: np.hypot(*(np.subtract(enemy(t), shot(t, angle)))), bounds=(LAUNCH_TIME, 250),
                          method="bounded", options={"xatol": 1e-9})
    return res.fun


p = argparse.ArgumentParser()
p.add_argument("--angles", type=float, nargs=3, default=[20.0, 30.0, 40.0])
p.add_argument("--solve", action="store_true")
a = p.parse_args()

angles = a.angles
if a.solve:
    best = minimize_scalar(miss_distance, bounds=(0, 89), method="bounded", options={"xatol": 1e-10}).x
    print(f"intercept angle {best:.6f} deg (miss distance {miss_distance(best):.3g} m)")
    angles = [float(best)] * 3
hit, t_end, closest = play(angles)
print(f"angles {angles}: {'intercepted' if hit else 'Fail'} at t = {t_end:.2f} s; closest approaches {closest} m")

tt = np.linspace(0, t_end, 500)
plt.plot(*enemy(tt), "o", markersize=2, label="enemy")
for k, ang in enumerate(angles):
    ts = tt[tt >= LAUNCH_TIME]
    plt.plot(*shot(ts, ang), "*", markersize=2, label=f"shot {k + 1} ({ang:.3f} deg)")
plt.plot(0, 0, "d", CITY, 0, "p")
plt.axis([-1, 3.5e5, -1, 1e5])
plt.legend()
plt.show()
