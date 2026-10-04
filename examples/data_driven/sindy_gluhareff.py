"""SINDy model of the Gluhareff pressure-jet engine pressures (``SINDy_Gluhareff_v_0``).

Needs ``Gluhareff_massAveDat_adiab.mat`` (variables ``P1, P2, P3, Pexit,
TimeMilisec``), which is not in the repository::

    python sindy_gluhareff.py Gluhareff_massAveDat_adiab.mat --second-order

Data are split 70/15/15 into training, validation and test; training
starts at sample 400 (steady operation). ``--second-order`` augments the
state with its TV-regularised derivatives, as in the final MATLAB version.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import SINDy
from unicodes.io import load_mat
from unicodes.numerics import tv_derivative

p = argparse.ArgumentParser()
p.add_argument("data")
p.add_argument("--second-order", action="store_true")
p.add_argument("--threshold", type=float, default=0.5)
p.add_argument("--ridge", type=float, default=0.1)
a = p.parse_args()

d = load_mat(a.data)
state = np.column_stack([np.ravel(d[k]) for k in ("P1", "P2", "P3", "Pexit")])
names = ["P1", "P2", "P3", "Pexit"]
t_ms = np.ravel(d["TimeMilisec"])
dt = (t_ms[1] - t_ms[0]) / 1000


def deriv(x):
    return tv_derivative(x, 10, 2e-5, dt, ep=1e12, scale="large")


if a.second_order:
    state = np.column_stack([state, np.column_stack([deriv(s) for s in state.T])])
    names += [n + "_dot" for n in names]
n = len(state)
i_train, i_val = int(0.7 * n), int(0.7 * n) + int(0.15 * n)
train = state[400:i_train]
dX = np.column_stack([deriv(s) for s in train.T])
model = SINDy(("const", "linear", "sin"), threshold=a.threshold, ridge=a.ridge).fit(train, dX)
print("\n".join(model.equations(names)))

for label, sl in (("training", slice(400, i_train)), ("validation", slice(i_train, i_val))):
    t = (t_ms[sl] - t_ms[sl][0]) / 1000
    sim = model.simulate(state[sl][0], t, method="LSODA")
    fig, axes = plt.subplots(2, 1, sharex=True)
    axes[0].plot(t_ms[sl][: len(sim)], sim)
    axes[0].set_title(f"model on {label} region")
    axes[1].plot(t_ms[sl], state[sl])
    axes[1].set_title(f"measured {label} data")
plt.show()
