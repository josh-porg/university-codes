"""AE 546 Lab 3: flat-plate boundary-layer thickness from a pitot rake.

Expects ``Lab3_data(1).xlsx``: first column labels, then per experiment a
static head, total head and 10 rake heads (cm of manometer fluid). The rake
pitots are 0.188333 cm apart, 0.3 m from the leading edge.
"""

import sys

import numpy as np

from unicodes.aero import boundary_layer as bl
from unicodes.aero.pressure import manometer_velocity

RHO_AIR, RHO_MANO, MU = 1.225, 826.0, 1.78e-5
X, DY = 0.3, 0.188333e-2


def analyse(h_s, h_t, h):
    v = manometer_velocity(h * 0.01, h_s * 0.01, RHO_MANO, RHO_AIR)
    vf = manometer_velocity(h_t * 0.01, h_s * 0.01, RHO_MANO, RHO_AIR)
    vf = np.where(vf > 0.99 * v.max(axis=0), v.max(axis=0), vf)
    Re = RHO_AIR * vf * X / MU
    ratio = v / vf
    edge = np.argmax(ratio >= 0.99, axis=0)  # first pitot at 99% of free stream
    delta_e = np.where(ratio.max(axis=0) >= 0.99, DY * (edge + 1), 10 * DY)
    return dict(
        V=vf,
        Re=Re,
        delta_laminar=bl.laminar_thickness(X, Re),
        delta_turbulent=bl.turbulent_thickness(X, Re),
        delta_experiment=delta_e,
        x_transition_start=bl.transition_distance(3e5, vf, MU / RHO_AIR),
        x_transition_end=bl.transition_distance(1e6, vf, MU / RHO_AIR),
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python lab3_boundary_layer.py 'Lab3_data(1).xlsx'")
    import pandas as pd

    d = pd.read_excel(sys.argv[1], header=None).to_numpy()[:, 1:].astype(float)
    for k, v in analyse(d[0], d[1], d[2:]).items():
        print(f"{k:20s}", np.array2string(v, precision=4))
