"""MATH 526: DMD of hospital-diagnosis data (``MedicalDataAnalysis/DataLoader``, ``DiagnosisDMD``).

Reads the MIMIC-IV ``diagnoses_icd.csv`` (columns ``subject_id, hadm_id,
seq_num, icd_code, icd_version``; not in the repository), one-hot encodes
each admission's ICD codes, pairs every patient's consecutive admissions
(``X`` = earlier visit, ``X'`` = next visit) and runs DMD on the pairs::

    python diagnosis_dmd.py diagnoses_icd.csv --threshold 4

The ten largest entries of the leading DMD mode name the diagnoses that
most strongly carry over between visits. (The MATLAB looked the codes up
in column 4 of the csv after building the mapping from column 3; the same
code column is used for both here.)
"""

import argparse
import csv
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

from unicodes.decomposition import exact_dmd

p = argparse.ArgumentParser()
p.add_argument("csv")
p.add_argument("--threshold", type=float, default=4.0, help="singular-value truncation threshold")
a = p.parse_args()

visits = defaultdict(set)
with open(a.csv, newline="") as f:
    for row in csv.DictReader(f):
        visits[(int(row["subject_id"]), int(row["hadm_id"]))].add(row["icd_code"].strip())
codes = sorted(set().union(*visits.values()))
index = {c: i for i, c in enumerate(codes)}
by_patient = defaultdict(list)
for (subject, hadm) in sorted(visits):
    by_patient[subject].append(hadm)

pairs = [(s, h0, h1) for s, adm in by_patient.items() for h0, h1 in zip(adm[:-1], adm[1:])]
X = np.zeros((len(codes), len(pairs)))
Xp = np.zeros_like(X)
for k, (s, h0, h1) in enumerate(pairs):
    X[[index[c] for c in visits[(s, h0)]], k] = 1
    Xp[[index[c] for c in visits[(s, h1)]], k] = 1
print(f"{len(visits)} admissions, {len(codes)} distinct codes, {len(pairs)} repeat-visit pairs")

lam, modes, S = exact_dmd(X, Xp, threshold=a.threshold)
lead = modes[:, np.argmax(np.abs(lam))]
lead = lead.real * np.sign(lead.real[np.argmax(np.abs(lead.real))])
top = np.argsort(lead)[::-1][:10]
print("leading eigenvalue", lam[np.argmax(np.abs(lam))])
print("key diagnoses:", [codes[i] for i in top])

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
axes[0].plot(S, "x")
axes[0].set(title="singular values", yscale="log")
axes[1].plot(lam.real, lam.imag, "x")
axes[1].set(title="DMD eigenvalues", xlabel="Re", ylabel="Im")
axes[2].plot(np.sort(lead)[::-1], "x")
axes[2].set(title="leading mode entries (sorted)")
plt.show()
