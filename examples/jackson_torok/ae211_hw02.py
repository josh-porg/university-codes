"""Jackson Torok, AE 211 homework 2 (``Torok_Jackson_HW2_MATLAB``): projectile drop, point-to-line
distance, rounding, matrix assembly (4.1-4.11).

Problem 4.5 reads ``ace_data.dat`` (year, ACE, tropical storms, hurricanes,
major hurricanes), which is not in the repository; pass it with
``--ace-data``. Fix: the point-to-line distance in 1b used ``zA + z0``
instead of ``zA - z0``.
"""

import argparse

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--ace-data", help="ace_data.dat for problem 4.5")
a = p.parse_args()

# 1a
t = np.arange(0, 5.01, 0.5)
print("1a: y =", 100 - 0.5 * 32.2 * t**2)
# 1b: distance from A to the line through P0 with direction (a, b, c)
d_vec = np.array([0.6, 0.5, 0.7])
P0, A = np.array([-4, -2, -3]), np.array([2, -3, 1])
r = A - P0
dA0 = np.linalg.norm(r)
angle = np.arccos(r @ d_vec / (dA0 * np.linalg.norm(d_vec)))
print(f"1b: distance = {dA0 * np.sin(angle):.6f} (check {np.linalg.norm(np.cross(r, d_vec)) / np.linalg.norm(d_vec):.6f})")
# 1c (MATLAB round(x, -3) rounds to thousands)
orig = 316501.673
print(f"1c: {round(orig, 2)}, {round(orig, -3):.0f}")
# 4.1
amat = np.array([[15, 3, 22], [3, 8, 5], [14, 3, 82]])
bmat = np.array([[1], [5], [6]])
cmat = np.array([[12, 18, 5, 2]])
dmat = amat[:, [2]]
print("4.1: e =\n", np.hstack([bmat, dmat]), "\n f =", np.vstack([bmat, dmat]).ravel(),
      "\n g =\n", np.vstack([amat, cmat[:, :3]]), "\n h =", np.array([amat[0, 2], cmat[0, 1], bmat[1, 0]]))
# 4.3
times = np.arange(0, 25, 2)
tc = np.array([[84.3, 90, 86.7], [86.4, 89.5, 87.6], [85.2, 88.6, 88.3], [87.1, 88.9, 85.3], [83.5, 88.9, 80.3],
               [84.8, 90.4, 82.4], [85.0, 89.3, 83.4], [85.3, 89.5, 85.4], [85.3, 88.9, 86.3], [85.2, 89.1, 85.3],
               [82.3, 89.5, 89.0], [84.7, 89.4, 87.3], [83.6, 89.8, 87.2]])
print("4.3: max times", times[tc.argmax(0)], "min times", times[tc.argmin(0)])
# 4.5
if a.ace_data:
    ace = np.genfromtxt(a.ace_data)
    ace = ace[~np.isnan(ace).any(axis=1)]  # drop a header line
    cols = ace[:, 1:5]
    print("4.5: years of the maxima", ace[cols.argmax(0), 0], "\n means", cols.mean(0), "\n medians", np.median(cols, 0))
    print(" sorted by ACE:\n", ace[np.argsort(-ace[:, 1])])
else:
    print("4.5: skipped (pass --ace-data ace_data.dat)")
# 4.8
pressure = np.arange(0, 100001, 10000)[:, None]
print("4.8: column heights [mercury, water]\n", pressure / (np.array([13560, 1000]) * 9.81))
# 4.10
print("4.10: zero matrices of shapes", np.zeros_like(amat).shape, np.zeros_like(bmat).shape, np.zeros_like(cmat).shape)
# 4.11 (MATLAB magic(6))
mag = np.array([[35, 1, 6, 26, 19, 24], [3, 32, 7, 21, 23, 25], [31, 9, 2, 22, 27, 20],
                [8, 28, 33, 17, 10, 15], [30, 5, 34, 12, 14, 16], [4, 36, 29, 13, 18, 11]])
print("4.11: row sums", mag.sum(1), "column sums", mag.sum(0), "diagonals", [np.trace(mag), np.trace(np.flipud(mag))])
