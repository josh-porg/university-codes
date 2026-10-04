"""AE 445 Homework 14/15, problem 4.4: drag polar of a wing from (alpha, C_L, C_D) data
(``HW_14_and_15_problem_4_4``).

Pass ``Prob_4_4_Data.xlsx`` (or a CSV with the same three columns). The
MATLAB only plotted C_D against C_L and C_L^2; the straight-line fit of
``C_D = C_D0 + C_L^2 / (pi A e)`` is added to read off C_D0 and e for
A = 10.32.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("data")
p.add_argument("--aspect-ratio", type=float, default=10.32)
a = p.parse_args()

if a.data.endswith((".xlsx", ".xls")):
    import pandas as pd  # reading .xlsx needs pandas + openpyxl

    data = pd.read_excel(a.data).to_numpy(float)
else:
    data = np.loadtxt(a.data, delimiter=",", skiprows=1)
alpha, cl, cd = data[:, 0], data[:, 1], data[:, 2]
k, cd0 = np.polyfit(cl**2, cd, 1)
print(f"C_D0 = {cd0:.5f}, k = {k:.5f}, e = {1 / (np.pi * a.aspect_ratio * k):.3f}")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(cl, cd, "o")
ax1.set(xlabel="C_L", ylabel="C_D")
ax2.plot(cl**2, cd, "o", cl**2, cd0 + k * cl**2, "-")
ax2.set(xlabel="C_L^2", ylabel="C_D")
plt.show()
