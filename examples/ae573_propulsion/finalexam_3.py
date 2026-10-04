"""AE 573 final exam problem 3 (``FinalExam_3``).

turbofan

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'M1 = v1/a1',
    'M1 = Ccz1/(gamma*R*T)^.5',
    'U = w*r',
    'Rdeg = 1- (ct1-ct2)/2*U',
    'Ct2 = 2*U*(1-Rdeg)',
    'psi = ct2/U-ct1/U',
    'Tt1/T1 = 1+(r-1)/2*M1^2',
    'Tt2/Tt1 = 1 +U*(ct2-ct1)/(cp*Tt1)',
    'T2 = Tt2-c2^2/(2*cp)',
    'a2 = (gamma*R*T2)',
    'M2=c2/a2',
    'c3 = c2',
    'Tt3 = Tt2',
    'T3=Tt3-c3/2*cp',
    'a3 = (gamma*R*T3)^.5',
    'M3 = c3/a3',
    'cp = cv + R',
    'cp = gamma * R / (gamma - 1)',
    'cv = R / (gamma - 1)',
    'v_2 = c2',
    'v_1 = cz1',
    'M2= M_2',
    'Tt2 =Tt_2',
    'pt_1 = pt1',
    'Tt1 =Tt_1',
    'pt_2 = pt2',
]

# (templates, token, values): each template is repeated with ``token`` -> ``_<value>``
REPEATED = [
    (
        [
            'M_st = V_st / a_st',
            'a_st = (gamma*R*T_st)^.5',
            'pt_st / p_st = (1 + (gamma-1)/2 * M_st^2)^(gamma/(gamma-1))',
            'rhot_st / rho_st = (1 + (gamma-1)/2 * M_st^2)^(1/(gamma-1))',
            'Tt_st / T_st = (1 + (gamma-1)/2 * M_st^2)^1',
            'p_st = rho_st * R * T_st',
        ],
        '_st',
        [0, 1, 2, 3, 13, 19],
    ),
]

K = {}
K["r"] = .4
K["cz1"] = 150
K["sigmar"] = 1.5
K["w"] = 774.926
K["Rdeg"] = .75
K["wtilde"] = .05
K["gamma"] = 1.4
K["gamma_1"] = 1.4
K["gamma_2"] = 1.4
K["R"] = 287
K["Pt1"] = 100e3
K["cx1"] = 150
K["c3"] = 150
K["Tt_1"] = 287

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
