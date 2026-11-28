"""R. Bramlette (University of Kansas, 2010), ``Bramlette_ku_0099D_14625_DATA_3gn``: sizing a
Gluhareff pressure-jet engine from its injector.

Propane at 74.7 psia and 420 deg F is injected through a 0.040 in.
orifice (efficiency 0.9097). The choked injector flow and the kinetic
energy of the jet set the three entrainment stages: stage k carries
``phi_k`` times the stoichiometric air/fuel ratio (15.5) of air per unit of
fuel and receives ``eta_i^k`` of the jet kinetic energy (inlet efficiency
0.6747). Rule-of-mixtures density, gas constant and gamma give each stage's
area, sound speed and (quarter- or half-wave) length at the operating
frequency, which scales from the G8-2-15 engine (639 Hz at a 0.125 in.
injector) with the injector diameter. Trend ratios then size the
combustor and nozzle. English units, as in the original.

Fix: the injector exit density used the feed pressure in psia with the
expanded exit temperature (``rho0 = Pf / (R T0)``); the jet is expanded to
ambient pressure, so ``rho0 = 144 Pamb / (R T0)`` in slug/ft^3 here (28
times the MATLAB value). Through the rule of mixtures that raises the
stage densities by 8-22 % and shrinks the inlet diameters from 0.320,
0.591 and 0.729 in. to 0.290, 0.565 and 0.701 in. (stage 3 length is
unchanged; stage 1 and 2 lengths shift by 2-5 %). The injector column shows
the exit state (ambient pressure). The stage static pressures and first
stage temperature carry the author's ``***CHECK THIS***`` notes and are
kept as written.
"""

import numpy as np

# design inputs
ID0, Pf, Tf_F, Tflame_F = 0.0400, 74.700, 420.00, 1670.0
N0, Ni = 0.9097, 0.6747
phi = {1: 0.3807, 2: 0.8068, 3: 0.9429}
AF_stoich = 15.500
F15, ID15 = 639.0, 0.125
# operating conditions
Pamb, Tamb_F, g = 14.700, 70.000, 32.178
y_prop, y_air, R_prop, R_air = 1.15, 1.40, 1130.0, 1716.0

A0 = np.pi * (ID0 / 2) ** 2  # in^2
Tf, Tamb, Tflame = Tf_F + 459.67, Tamb_F + 459.67, Tflame_F + 459.67
rho_amb = Pamb * 144 / (R_air * Tamb)

# injector
mdot0 = N0 * Pf * A0 * np.sqrt(y_prop / (R_prop * Tf) * (2 / (y_prop + 1)) ** ((y_prop + 1) / (y_prop - 1)))  # slug/s
T0 = Tf * (Pamb / Pf) ** ((y_prop - 1) / y_prop)
a0 = np.sqrt(y_prop * R_prop * T0)
rho0 = Pamb * 144 / (R_prop * T0)
V0 = N0 * a0
Q0 = 0.5 * rho0 * V0**2 / 144
KE0 = 0.5 * mdot0 * V0**2
Fx = F15 * ID15 / ID0
stages = {0: dict(eta=N0, mdot=mdot0, rho=rho0, L=0.0, A=A0, ID=ID0, V=V0, gamma=y_prop, R=R_prop, T=T0, a=a0,
                  P=Pamb, Q=Q0, f=1.0, KE=KE0)}

# inlets
for k in (1, 2, 3):
    f = phi[k] * AF_stoich
    mdot = f * mdot0
    eta = Ni**k
    KE = eta * KE0
    V = np.sqrt(2 * KE / mdot)
    rho = rho_amb * (1 - 1 / f) + rho0 / f
    A = 144 * mdot / (rho * V)
    Q = 0.5 * rho * V**2 / 144
    P = Pamb - Q
    R = R_air * (1 - 1 / f) + R_prop / f
    gam = y_air * (1 - 1 / f) + y_prop / f
    if k == 1:
        T = Tamb * (1 - 1 / f) + (P * 144 / (rho * R)) / f
    elif k == 2:
        T = P * 144 / (rho * R)
    else:
        T = Tflame  # assumed
    a = np.sqrt(gam * R * T)
    L = a / ((2 if k < 3 else 4) * Fx)
    stages[k] = dict(eta=eta, mdot=mdot, rho=rho, L=L, A=A, ID=2 * np.sqrt(A / np.pi), V=V, gamma=gam, R=R, T=T, a=a,
                     P=P, Q=Q, f=f, KE=KE)

# trend-based sizing
L2, L3 = stages[2]["L"], stages[3]["L"]
Ln = 2 * L2
IDcc = L3 / 2
Lcc = (3.6988 * ID0 + 0.8307) * L2 * 1.15
IDn = 2 / 3 * IDcc
ID3f = 1.25 * stages[3]["ID"]

print("The engine has been sized based on:")
print(f"Injector inner diameter = {ID0:.4f} in.\nInjector feed press.    = {Pf:.1f} psi\n"
      f"Injector feed temp.     = {Tf_F:.1f} deg F\nInjector efficiency     = {100 * N0:.2f} %\n"
      f"Inlet efficiency        = {100 * Ni:.2f} %\nEquivalency ratios      = {phi[1]:.4f} {phi[2]:.4f} {phi[3]:.4f}\n"
      f"Stoich. air/fuel ratio  = {AF_stoich:.2f}\n" + "=" * 45)
print(f"\nOperating frequency = {Fx:.1f} Hz\n")
for k, name in ((1, "First"), (2, "Second"), (3, "Third")):
    print(f"{name} stage inlet: inner diameter {stages[k]['ID']:.4f} in., length {12 * stages[k]['L']:.4f} in.")
print(f"Third stage flared ID {ID3f:.4f} in.")
print(f"Combustion chamber: inner diameter {12 * IDcc:.4f} in., length (with nose) {12 * Lcc:.4f} in.")
print(f"Nozzle: inner diameter {12 * IDn:.4f} in., length {12 * Ln:.4f} in.\n")
keys = ["eta", "mdot", "rho", "L", "A", "ID", "V", "gamma", "R", "T", "a", "P", "Q", "f", "KE"]
units = ["", "slug/s", "slug/ft3", "ft", "in2", "in", "ft/s", "", "ft lbf/slug R", "deg F", "ft/s", "psia", "psi", "", "ft lbf/s"]
print(f"{'':10s}" + "".join(f"{s:>13s}" for s in ("injector", "stage 1", "stage 2", "stage 3")))
for key, unit in zip(keys, units):
    vals = [stages[k][key] - (459.67 if key == "T" else 0) for k in range(4)]
    print(f"{key:6s}{unit:>14s}" + "".join(f"{v:13.5g}" for v in vals))
