"""Small stand-alone scripts from ``unsorted`` and elsewhere.

* ``Logistic``: bifurcation diagram of the logistic map obtained by iterating
  its inverse branch ``x = (b - sqrt(b (b - 4 x))) / (2 b)``.
* ``sumDidgetsOfPi``: sum of the first 10 000 digits of pi.
* ``weather_probability_advancing``: three-state Markov weather model.
* ``trackGearCalculator``: gear inches of a 48 x 15 track bike on 700c wheels.
* ``autumn 2023/orthogonal_vectors``: a unit vector orthogonal to two
  vectors (the MATLAB tried ``fsolve`` on symbolic equations; it is the
  normalised cross product).
"""

import matplotlib.pyplot as plt
import mpmath
import numpy as np

# Logistic map via the inverse iteration
betas = np.arange(0.001, 4.0, 0.002)
pts_b, pts_x = [], []
for b in betas:
    x = 0.5
    for _ in range(2000):  # discard the transient
        x = (b - np.sqrt(max(b * (b - 4 * x), 0))) / (2 * b)
    x_ss = x
    for _ in range(2000):
        x = (b - np.sqrt(max(b * (b - 4 * x), 0))) / (2 * b)
        pts_b.append(b)
        pts_x.append(x)
        if abs(x - x_ss) < 1e-3:
            break
fig, ax = plt.subplots(facecolor="k")
ax.plot(pts_x, 4 - np.array(pts_b), ".", ms=0.8, color=(1, 0.1, 1))
ax.set(facecolor="k", xlabel="x", ylabel="beta (4 at top)", ylim=(0, 1))

# Sum of the digits of pi
mpmath.mp.dps = 10010
digits = mpmath.nstr(mpmath.pi, 10001, strip_zeros=False).replace(".", "")[:10001]
print("sum of the first 10001 digits of pi (including the leading 3):", sum(map(int, digits)))

# Markov weather model
T = np.array([[0.5, 0.5, 0.25], [0.25, 0, 0.25], [0.25, 0.5, 0.5]])
w, v = np.linalg.eig(T)
steady = np.real(v[:, np.argmin(abs(w - 1))])
print("steady-state weather probabilities:", steady / steady.sum(), "(expected [0.4, 0.2, 0.4])")
state = np.array([0.0, 0, 1])
hist = []
for _ in range(100):
    hist.append(state)
    state = T @ state
fig, ax = plt.subplots()
ax.plot(np.array(hist))
ax.set(xlabel="day", ylabel="probability", title="Markov weather model")

# Track gear
diam_in = 700 / 25.4
print(f"gear ratio {48 / 15:.3f}, gear inches {48 / 15 * diam_in:.1f}, wheel circumference {np.pi * diam_in:.1f} in")

# Orthogonal unit vector
u, v_ = np.array([1.0, 2, 3]), np.array([-1.0, 0, 2])
w_ = np.cross(u, v_)
w_ /= np.linalg.norm(w_)
print("unit vector orthogonal to", u, "and", v_, ":", w_, "dots", u @ w_, v_ @ w_)
plt.show()
