"""AE 573 homework 8 part 3 (``H_8_part_3``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'Tt_1 = Tt_0',
    'p_9 = p_0',
    'p_19 = p_0',
    'alpha = mdot_fan / mdot_core',
    'tau_c = pi_c^ ((gamma-1)/(gamma*e_c))',
    'tau_f = pi_f^ ((gamma-1)/(gamma*e_f))',
    'P_fan = mdot_fan *cp * (Tt_1 - Tt_13)',
    'P_compressor = mdot_core *cp * (Tt_2 - Tt_3)',
    'mdot_core = mdot_2',
    'mdot_fan = mdot_2',
    'pi_d = pt_1/pt_0',
    'pi_f = pt_13/pt_1',
    'pi_f = pt_2/pt_1',
    'pi_fn = pt_19/pt_13',
    'pi_c = pt_3/pt_2',
    'tau_c = Tt_3 / Tt_2',
    'tau_f = Tt_2 / Tt_1',
    'tau_f = Tt_13 / Tt_1',
    'cp = cv + R',
    'cp = gamma * R / (gamma - 1)',
    'cv = R / (gamma - 1)',
]

# (templates, token, values): each template is repeated with ``token`` -> ``_<value>``
REPEATED = [
    (
        [
            'tau_x = pi_x^ ((gamma -1)/(gamma *e_x ))',
        ],
        '_x',
        ['d', 'f', 'c', 'fn', 'n'],
    ),
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
K["pi_c"] = 32
K["pi_d"] = .95
K["pi_f"] = 1.5
K["pi_fn"] = .98
K["e_c"] = .9
K["e_f"] = .9
K["M_0"] = .8
K["p_0"] = 30E3
K["T_0"] = 228
K["gamma"] = 1.4
K["R"] = 287
K["mdot_fan"] = 100
K["mdot_core"] = 20

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
