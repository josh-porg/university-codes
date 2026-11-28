"""AE 360 high-altitude balloon payload data (``ae_360_analysis``, first part).

Reads ``Flight Data.xlsx`` (time, lat, long, speed, course, alt, temp, press,
photo, RH), applies the sensor conversions from the MATLAB and plots each
channel against time and course.
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_excel(sys.argv[1]).iloc[1:]  # first row was a test point
t = (pd.to_datetime(df.iloc[:, 0]) - pd.to_datetime(df.iloc[0, 0])).dt.total_seconds().to_numpy()
names = ["lat", "long", "speed", "course", "alt", "temp", "press", "photo", "RH"]
data = df.iloc[:, 1:10].to_numpy(dtype=float)
data[:, 5] = data[:, 5] / 256 * 500 + 273  # temperature, sensor counts -> K (the MATLAB applied this to a row, not the column)
data[:, 6] = data[:, 6] * 0.8665 - 67.73  # pressure calibration guessed in the MATLAB
data[:, 8] = (data[:, 8] * 5 / 256 - 0.8) / 0.031  # relative humidity %
fig, axes = plt.subplots(3, 3, figsize=(12, 8))
for ax, name, col in zip(axes.flat, names, data.T):
    ax.plot(t, col)
    ax.set_title(name)
plt.figure()
ax = plt.subplot(projection="polar")
ax.plot(np.deg2rad(data[:, 3]), data[:, 7], ".")
ax.set_title("course vs photo sensor")
plt.show()
