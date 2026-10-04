"""Jackson Torok, AE 211 homework 7 (``Torok_Jackson_HW7_MATLAB``): logical functions and selection
(8.3-8.20): engine test screening, piecewise functions, sea-ice extent, menus and parking fees.

The MATLAB menus and ``input`` prompts are command-line options here.
Fixes: the ice-extent averages were drawn as two crossed lines
(``plot([1979,1979;2010,2010],[15.52,6.51;15.52,6.51])``) and labelled the
wrong way round; they are horizontal lines at the March (15.52) and
September (6.51) means. Parking rates as in the MATLAB: long term $2 for the
first hour plus $1 per further hour (at most $9 a day), $9 a day (at most
$60 a week), $60 a week; short term $2 for the first 30 minutes plus $1 per
further 20 minutes (at most $32 a day).
"""

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--xy", type=float, nargs=2, default=[3, 2], help="8.12 x and y")
p.add_argument("--arcsin", type=float, default=0.5, help="8.13 argument")
p.add_argument("--major", type=int, choices=range(1, 6), default=5, help="8.18 menu choice")
p.add_argument("--star", type=int, choices=[1, 2], default=1, help="8.19: five or six points")
p.add_argument("--lot", choices=["long", "short"], default="short", help="8.20 parking lot")
p.add_argument("--stay", type=float, nargs=3, default=[0, 1, 3],
               help="8.20 long: weeks days hours; short: days hours minutes")
a = p.parse_args()

# 8.3
data = np.array([[1530, 116, 45, 110], [1240, 114, 42, 115], [2380, 118, 41, 120], [1470, 124, 38, 95],
                 [3590, 126, 61, 118]])
temp_ok = (data[:, 1] > 115) & (data[:, 1] < 125)
hum_ok = (data[:, 2] > 40) & (data[:, 2] < 60)
pres_ok = (data[:, 3] > 100) & (data[:, 3] < 200)
for label, ok in (("Temp", temp_ok), ("Humidity", hum_ok), ("Pressure", pres_ok), ("Every Test", temp_ok & hum_ok & pres_ok)):
    print(f"Engines {data[ok, 0]} passed {label}")
print(f"The total pass rate was {100 * np.mean(temp_ok & hum_ok & pres_ok):.2f}%")


# 8.6
def g(x):
    return np.where((x >= -np.pi) & (x <= np.pi), np.cos(x), -1.0)


x = np.linspace(-2 * np.pi, 2 * np.pi, 500)
plt.figure(num="Question 8.6")
plt.plot(x, g(x))
# 8.9
ice = np.array([[1979, 7.19, 16.48], [1980, 7.83, 16.15], [1981, 7.24, 15.65], [1982, 7.44, 16.17],
                [1983, 7.51, 16.13], [1984, 7.10, 15.65], [1985, 6.91, 16.09], [1986, 7.53, 16.10],
                [1987, 7.47, 15.99], [1988, 7.48, 16.16], [1989, 7.03, 15.54], [1990, 6.23, 15.90],
                [1991, 6.54, 15.52], [1992, 7.54, 15.50], [1993, 6.50, 15.90], [1994, 7.18, 15.62],
                [1995, 6.12, 15.35], [1996, 7.87, 15.16], [1997, 6.73, 15.61], [1998, 6.55, 15.69],
                [1999, 6.23, 15.45], [2000, 6.31, 15.30], [2001, 6.74, 15.64], [2002, 5.95, 15.46],
                [2003, 6.3, 15.52], [2004, 6.04, 15.08], [2005, 5.56, 14.77], [2006, 5.91, 14.45],
                [2007, 4.29, 14.66], [2008, 4.72, 15.27], [2009, 5.38, 15.16], [2010, 4.92, 15.14]])
plt.figure(num="Question 8.9")
plt.plot(ice[:, 0], ice[:, 1], label="September")
plt.plot(ice[:, 0], ice[:, 2], label="March")
plt.hlines([15.52, 6.51], 1979, 2010, linestyles="--", colors=["C1", "C0"], label=None)
plt.legend(loc="lower right")
plt.axis([1979, 2010, 0, 20])
print("8.9: September above average in", ice[ice[:, 1] > 6.51, 0].astype(int))
print("     March above average in", ice[ice[:, 2] > 15.52, 0].astype(int))
# 8.12
x_, y_ = a.xy
print("8.12:", "x>y" if x_ > y_ else "x<=y")
# 8.13
print("8.13:", math.asin(a.arcsin) if -1 < a.arcsin < 1 else "Arcsin is non-real at the current input")
# 8.18
credits = [130, 130, 122, 126.5, 129]
print(f"8.18: you'll need {credits[a.major - 1]:5.1f} credits to graduate")
# 8.19
fig = plt.figure(num="Question 8.19")
ax = fig.add_subplot(projection="polar")
if a.star == 1:
    theta = np.arange(np.pi / 2, 4.5 * np.pi + 1e-9, 4 * np.pi / 5)
    ax.plot(theta, np.ones_like(theta))
else:
    for start in (np.pi / 2, np.pi / 6):
        theta = np.arange(start, 2 * np.pi + start + 1e-9, 2 * np.pi / 3)
        ax.plot(theta, np.ones_like(theta))
# 8.20
if a.lot == "long":
    w, d, h = a.stay
    total = w * 60 + min(d * 9, 60) + min(2 + (h - 1), 9)
else:
    d, h, m = a.stay
    minutes = h * 60 + m
    total = d * 32 + min(max(2 + math.ceil((minutes - 30) / 20), 2), 32)
print(f"8.20: parking fee ${total:g}")
plt.show()
