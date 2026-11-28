"""Jackson Torok, AE 211 homework 4 (``Torok_Jackson_HW4_MATLAB``): plotting (5.14-5.31) — projectile
paths, polar plots, Arrhenius rates, grade histograms and pie charts, 3-D curves and surfaces.

The 5.14 legend labelled all three launch angles pi/2; here they are
labelled pi/2, pi/4 and pi/6.
"""

import matplotlib.pyplot as plt
import numpy as np

# 5.14
theta = np.array([np.pi / 2, np.pi / 4, np.pi / 6])[:, None]
t = np.arange(0, 20.01, 0.1)
x = t * 100 * np.cos(theta)
y = t * 100 * np.sin(theta) - 0.5 * 9.8 * t**2
fig, ax = plt.subplots(2, 2, num="Question 5.14")
ax[0, 0].plot(t, x[1])
ax[0, 0].set(title="Horizontal distance vs time", xlabel="Time (s)", ylabel="m")
ax[0, 1].plot(t, y[1])
ax[0, 1].set(title="Vertical distance vs time", xlabel="Time (s)", ylabel="m")
ax[1, 0].plot(x[1], y[1])
ax[1, 0].set(title="Vertical vs horizontal distance")
for k, (style, label) in enumerate((("--b", r"$\theta=\pi/2$"), ("-g", r"$\theta=\pi/4$"), (":k", r"$\theta=\pi/6$"))):
    ax[1, 1].plot(x[k], y[k], style, label=label)
ax[1, 1].legend()
# 5.15
th = np.linspace(0, 2 * np.pi, 101)
fig = plt.figure(num="Question 5.15")
for k, r in enumerate((np.sin(th) ** 2 + np.cos(th) ** 2, np.sin(th), np.exp(th / 5), np.sinh(th))):
    fig.add_subplot(2, 2, k + 1, projection="polar").plot(th, r)
# 5.20
T = np.arange(300, 1001)
k_rate = 10 * np.exp(-1000 / (8.314 * T))
fig, ax = plt.subplots(2, 1, num="Question 5.20")
ax[0].plot(T, k_rate)
ax[0].set(title="Reaction rate k vs temperature (K)")
ax[1].plot(1 / T, np.log(k_rate))
ax[1].set(title="log(k) vs 1/T")
# 5.21 / 5.22
G = np.array([68, 83, 61, 70, 75, 82, 57, 5, 76, 85, 62, 71, 96, 78, 76, 68, 72, 75, 83, 93])
edges = [0, 60, 70, 80, 90, 100]
fig, ax = plt.subplots(2, 2, num="Question 5.21")
ax[0, 0].bar(range(G.size), np.sort(G))
ax[0, 0].set(title="Part a")
ax[0, 1].hist(G)
ax[0, 1].set(title="Part b")
ax[1, 0].hist(G, bins=edges)
ax[1, 0].set(title="Part c")
ax[1, 1].hist(G, bins=edges, density=True)
ax[1, 1].set(title="Part d")
n, _ = np.histogram(G, edges)
print("5.22: grade counts E-A", n)
plt.figure(num="Question 5.22")
plt.pie(n, labels=list("EDCBA"))
# 5.26
fig, ax1 = plt.subplots(num="Question 5.26")
ax1.plot(t, x[1], "b")
ax1.set(ylabel="Distance (m)", xlabel="Time (s)")
ax2 = ax1.twinx()
ax2.plot(t, y[1], "r")
ax2.set(ylabel="Height (m)")
# 5.29
x1 = np.arange(0, 20 * np.pi + 1e-9, np.pi / 100)
fig = plt.figure(num="Question 5.29")
fig.add_subplot(2, 2, 1).plot(x1, x1 * np.sin(x1))
fig.add_subplot(2, 2, 2, projection="3d").plot(x1, x1 * np.sin(x1), x1 * np.cos(x1))
fig.add_subplot(2, 1, 2, projection="polar").plot(x1, x1 * np.sin(x1))
# 5.31
X, Y = np.meshgrid(np.arange(-5, 5.01, 0.5), np.arange(-5, 5.01, 0.5))
Z = np.sin(np.sqrt(X**2 + Y**2))
fig = plt.figure(num="Question 5.31")
fig.add_subplot(2, 3, 1, projection="3d").plot_wireframe(X, Y, Z)
fig.add_subplot(2, 3, 2, projection="3d").plot_surface(X, Y, Z)
fig.add_subplot(2, 3, 3, projection="3d").plot_surface(X, Y, Z, cmap="jet")
cs = fig.add_subplot(2, 3, 4).contour(X, Y, Z)
plt.clabel(cs)
ax3 = fig.add_subplot(2, 3, 5, projection="3d")
ax3.plot_surface(X, Y, Z, alpha=0.5)
ax3.contour(X, Y, Z, offset=-1)
plt.show()
