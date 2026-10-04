"""Electrical power budget of the sailplane by flight phase (``Electrical Systems/PowerDiagram``).

Reads the team's ``Load Tracker.csv`` (columns ElectricalLoadItem, Model,
ElectricalSystem, PreflightW, TakeoffClimbW, SoaringFlightW, LandingTaxiW;
not in the repository) and plots the power drawn by each item in each phase::

    python power_budget.py "Load Tracker.csv"
"""

import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("csv")
a = p.parse_args()

phases = ["PreflightW", "TakeoffClimbW", "SoaringFlightW", "LandingTaxiW"]
items, power = [], []
with open(a.csv, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row.get("ElectricalLoadItem"):
            items.append(row["ElectricalLoadItem"].strip())
            power.append([float(row[ph] or 0) for ph in phases])
power = np.array(power)
for ph, total in zip(phases, power.sum(axis=0)):
    print(f"{ph[:-1]:15s} {total:8.1f} W")
fig, ax = plt.subplots(figsize=(9, 5))
bottom = np.zeros(len(phases))
for name, row in zip(items, power):
    ax.bar([ph[:-1] for ph in phases], row, bottom=bottom, label=name)
    bottom += row
ax.set(ylabel="power (W)", title="Electrical load by flight phase")
ax.legend(fontsize=7, ncol=2)
plt.show()
