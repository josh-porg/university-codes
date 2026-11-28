"""Jackson Torok, AE 211 homework 9 (``Torok_Jackson_HW9_MATLAB``): data types (11.1-11.19) — single
versus double round-off, integer division, overflow, character codes, string
tables and multidimensional arrays.

Problem 11.17 loads ``test_results.mat`` (``test1year1`` ... ``test3year2``),
which is not in the repository; pass it with ``--test-results``.
"""

import argparse

import numpy as np

from unicodes.io import load_mat

p = argparse.ArgumentParser()
p.add_argument("--test-results")
a = p.parse_args()

# 11.1
n = np.arange(1, 10_000_001)
print(f"11.1: double sum {np.sum(1.0 / n):.15f}, single sum {np.cumsum(1 / n.astype(np.float32), dtype=np.float32)[-1]:.6f}"
      " (adding term by term in single precision stalls once 1/n drops below the round-off of the sum)")
# 11.2: MATLAB int8 arithmetic rounds 1 ./ int8(1:10) to 1, 1 (0.5 rounds up), 0, ...
nd = np.arange(1, 11)
int_recip = np.round(1 / nd + 1e-12).astype(np.int8)  # MATLAB rounds halves away from zero
print(f"11.2: expected {np.sum(1 / nd):.6f}; int8 terms {int_recip}, sum {int_recip.sum()}")
# 11.4
with np.errstate(over="ignore", invalid="ignore"):
    print("11.4: double (5+3i)^100 =", complex(5 + 3j) ** 100, "; single =", np.complex64(5 + 3j) ** 100,
          "(single overflows to inf first)")
# 11.6
print(f"11.6: '85' has {len('85')} characters; codes {ord('8')} and {ord('5')}")
# 11.7 / 11.8
names = ["Emeliann", "Isabel", "Kayleigh", "Michelle", "Matt"]
bd = np.array([[3, 6, 2002], [2, 18, 2002], [3, 22, 2002], [5, 20, 2002], [12, 4, 2001]])
for name, (m, d, y) in zip(names, bd):
    print(f"{name:>10s} {m:3d} {d:3d} {y:5d}")
# 11.15
ABC = np.stack([[[1, 2], [3, 4]], [[10, 20], [30, 40]], [[3, 6], [9, 12]]], axis=2)
print("11.15: column 1 of each page\n", ABC[:, 0, :], "\n row 2 of each page\n", ABC[1, :, :].T,
      "\n element (1, 2, 3):", ABC[0, 1, 2])
# 11.17
if a.test_results:
    m = load_mat(a.test_results)
    data = np.stack([np.stack([m[f"test{t}year{y}"] for t in (1, 2, 3)], axis=-1) for y in (1, 2)], axis=2)
    print("11.17: student 1, question 2, year 1, test 3:", data[0, 1, 0, 2])
    print(" student 1, question 1, test 2, both years:", data[0, 0, :, 1])
    print(" student 2, test 1, year 2:", data[1, :, 1, 0])
    print(" question 3, test 2:\n", data[:, 2, :, 1])
else:
    print("11.17: skipped (pass --test-results test_results.mat)")
# 11.19
print("11.19:", ["aluminum", "copper", "iron", "molybdenum", "cobalt"])
