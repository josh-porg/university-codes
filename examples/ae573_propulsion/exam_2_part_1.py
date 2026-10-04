"""AE 573 exam 2 part 1 (``Exam_2_part_1``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'NPR = pt_7/p_0',
    'pi_n = pt_9/pt_7',
    'NPR_crit = 1/pi_cn*((gamma+1)/2)^(gamma/(gamma-1))',
    'NPR_crit = pt_7/p_0',
    'M_9 = (2/(gamma-1)*(pt_9/p_9)^((gamma-1)/gamma)-1)',
    'p8rat = (1 + (gamma-1)/2 * M_8^2)^(gamma/(gamma-1))*((gamma+1)/2)^(gamma/(gamma-1))',
    'p8rat = p8/p0',
    'cp = cv + R',
    'cp = gamma * R / (gamma - 1)',
    'cv = R / (gamma - 1)',
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
        [7, 8, 9],
    ),
]

K = {}
K["NPR"] = 3.78
K["pi_cn"] = .98
K["Tt_9"] = 800
K["gamma"] = 1.33
K["R"] = 287
K["cp"] = 1156

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
