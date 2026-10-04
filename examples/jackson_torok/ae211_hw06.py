"""Jackson Torok, AE 211 homework 6 (``Torok_Jackson_HW6_MATLAB``): input and output (7.3-7.17).

The MATLAB asked for keyboard input (cone base and height, age, an array)
and picked two points on a plot with ``ginput``; these are command-line
options here, with the two points defaulting to opposite ends of the unit
circle.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--base-area", type=float, default=3.0)
p.add_argument("--height", type=float, default=4.0)
p.add_argument("--age", type=float, default=20)
p.add_argument("--array", type=float, nargs="+", default=[1, 2, 3, 4])
p.add_argument("--points", type=float, nargs=4, default=[1, 0, -1, 0], metavar=("X1", "Y1", "X2", "Y2"))
a = p.parse_args()

print(f"7.3: the volume of the cone is {a.base_area * a.height / 3:g}")
print(f"7.6: your age is {a.age:g}")
yen = np.arange(5, 126, 5)
print("7.12: yen -> dollars\n", np.column_stack([yen, yen * 0.0092])[:5], "...")
euro = np.arange(1, 61, 2)
print(" euro -> dollars\n", np.column_stack([euro, euro * 1.19])[:5], "...")
dol = np.arange(1, 11)
print(" dollar, euro, pound, yen\n", np.column_stack([dol, dol / 1.19, dol / 1.39, dol / 0.0092]))
print(f"7.7: you entered a {len(a.array)} column long array")
print("7.13:")
for row in zip(["Smith", "Jones", "Webb", "Anderson"], ["Fred", "Kathy", "Milton", "John"], [6, 22, 92, 45],
               [47, 66, 62, 72], [82, 140, 110, 190]):
    print("  {:10s} {:8s} {:4d} {:4d} {:4d}".format(*row))
# 7.17
t = np.arange(0, 2 * np.pi + 1e-9, np.pi / 100)
plt.plot(np.cos(t), np.sin(t))
x1, y1, x2, y2 = a.points
plt.plot([x1, x2], [y1, y2])
plt.axis("equal")
print(f"7.17: distance between the points = {np.hypot(x2 - x1, y2 - y1):.4f}")
plt.show()
