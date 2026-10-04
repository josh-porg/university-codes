"""Jackson Torok, AE 211 homework 8 (``Torok_Jackson_HW8_MATLAB``): repetition structures (9.5-9.18):
prime products, Fibonacci sequences and the golden ratio, an iterative square
root, the sine series, the Leibniz series for pi and a lake-level report.

The Fibonacci seeds and term counts were typed in at prompts; they are
options here. ``report.txt`` (months above the eight-year average) is
written next to the script unless ``--report`` says otherwise.
"""

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
from sympy import primerange

p = argparse.ArgumentParser()
p.add_argument("--fib-start", type=float, nargs=2, default=[1, 1])
p.add_argument("--terms", type=int, default=20)
p.add_argument("--report", default="report.txt")
a = p.parse_args()

# 9.5
primes = np.array(list(primerange(2, 101)))
print("9.5: products of consecutive primes", primes[:-1] * primes[1:])


# 9.6 / 9.7 (the for and while versions give the same sequence)
def fibonacci(start, count):
    f = list(start)
    while len(f) < count:
        f.append(f[-2] + f[-1])
    return np.array(f[:count])


f = fibonacci(a.fib_start, a.terms)
fig = plt.figure(num="Question 9.6/9.7")
fig.add_subplot(projection="polar").plot(np.arange(1, f.size + 1), f)
# 9.8
x = list(a.fib_start)
x.append(x[-2] + x[-1])
while abs(x[-1] / x[-2] - x[-2] / x[-3]) > 0.001:
    x.append(x[-2] + x[-1])
print(f"9.8: phi = {x[-1] / x[-2]:.6f}, last change {x[-1] / x[-2] - x[-2] / x[-3]:.2e}")


# 9.13
def my_sqrt(num, guess, tol):
    """Third-order iteration g <- g/8 (15 - y (10 - 3 y)) with y = g^2 / num."""
    g_old, y = guess, guess**2 / num
    g = g_old / 8 * (15 - y * (10 - 3 * y))
    while abs(g - g_old) >= tol:
        g_old, y = g, g**2 / num
        g = g_old / 8 * (15 - y * (10 - 3 * y))
    return g


print(f"9.13: my_sqrt(5) = {my_sqrt(5, 2, 1e-4):.8f}, error {abs(my_sqrt(5, 2, 1e-4) - math.sqrt(5)):.2e}")


# 9.15
def my_sin(x):
    s = [x, x - x**3 / math.factorial(3)]
    k = 2
    while abs(s[-1] - s[-2]) > 0.001:
        k += 1
        s.append(s[-1] + (-1) ** (k - 1) * x ** (2 * k - 1) / math.factorial(2 * k - 1))
    return s[-1]


print(f"9.15: my_sin(2) = {my_sin(2):.6f}, error {abs(my_sin(2) - math.sin(2)):.2e}")
# 9.16
pi_sum, n = 1.0, 1
while n < 3000:
    n += 1
    term = (-1) ** (n - 1) / (2 * n - 1)
    pi_sum += term
    if abs(term) < 0.001:
        break
print(f"9.16: Leibniz pi = {4 * pi_sum:.6f} after {n} terms, error {4 * pi_sum - math.pi:.2e}")
# 9.18
lake = np.array([
    [3590.66, 3614.17, 3622.14, 3620.55, 3636.91, 3604.42, 3578.69, 3593.57],
    [3590.66, 3612.05, 3620.16, 3614.95, 3635.28, 3601.47, 3575.55, 3592.23],
    [3589.77, 3610.43, 3619.41, 3610.73, 3635.33, 3598.96, 3574.76, 3591.02],
    [3594.09, 3611.26, 3620.50, 3611.93, 3635.76, 3596.53, 3577.56, 3590.18],
    [3610.81, 3629.09, 3625.96, 3623.13, 3636.83, 3599.44, 3589.38, 3597.27],
    [3631.05, 3640.49, 3638.82, 3648.98, 3633.90, 3600.07, 3609.19, 3613.54],
    [3633.00, 3641.14, 3636.52, 3660.86, 3628.45, 3594.17, 3608.05, 3612.62],
    [3629.55, 3637.50, 3634.55, 3655.34, 3623.62, 3589.64, 3605.82, 3609.07],
    [3626.90, 3635.37, 3633.66, 3653.01, 3621.56, 3591.25, 3605.53, 3606.01],
    [3623.82, 3633.52, 3634.08, 3650.27, 3619.46, 3590.88, 3605.57, 3606.44],
    [3621.90, 3631.10, 3630.31, 3645.67, 3615.10, 3587.90, 3601.87, 3605.47],
    [3617.89, 3626.22, 3626.54, 3639.75, 3609.82, 3584.43, 3597.75, 3600.80]])
years = np.arange(2008, 2016)
avg = lake.mean()
with open(a.report, "w") as fh:
    for col, year in enumerate(years):
        for month in np.flatnonzero(lake[:, col] > avg):
            fh.write(f"Month {month + 1:2d} in {year} went above the 8 year average\n")
print(f"9.18: average {avg:.2f}; months above it per year {(lake > avg).sum(0)} (listed in {a.report})")
plt.show()
