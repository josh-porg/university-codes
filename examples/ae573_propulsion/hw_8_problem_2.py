"""AE 573 homework 8 problem 2 (``HW_8_problem_2``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'TSFC = mdot_f / Fn',
    'TSFC = f / (Fn/mdot_0)',
    'eta_th = DeltaKEdot / pth',
    'pth = mdot_f * Qr',
    'DeltaKEdot = (mdot_0 + mdot_f) * V_9^2 / 2 - mdot_0 * V_0^2 /2',
    'eta_p = Fn * V_0 / DeltaKEdot',
    'eta_p = 2/(1+V_9/V_0)',
    '(1+f)*cp_4*Tt_4 - cp_3*Tt_3 = f * Qr',
    'Tt_1 = Tt_0',
    'p_9 = p_0',
    'p_19 = p_0',
    'p_13 = p_19',
    'Fn = (mdot_0 + mdot_f) * V_9 - mdot_0 * V_0 + (p_9-p_0)*A_9',
    'alpha = mdot_fan / mdot_core',
    'mdot_0 = mdot_core+mdot_fan',
    'tau_c = pi_c^ ((gamma-1)/(gamma*e_c))',
    'tau_f = pi_f^ ((gamma-1)/(gamma*e_f))',
    'P_fan = mdot_fan *cp * (Tt_1 - Tt_13)',
    'P_compressor = mdot_core *cp * (Tt_2 - Tt_3)',
    'mdot_core = mdot_3',
    'mdot_fan = mdot_13',
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
K["V_0"] = 0
K["p_0"] = 101325
K["A_9"] = 1
K["pt_2"] = 14.7 * 6.895*10**3
K["pt_13"] = 26 * 6.895*10**3
K["pt_2.5"] = 63 * 6.895*10**3
K["pt_3"] = 200 * 6.895*10**3
K["pt_4"] = 190 * 6.895*10**3
K["pt_5"] = 28 * 6.895*10**3
K["Tt_2"] = 59 + 460
K["Tt_13"] = 170 + 460
K["Tt_2.5"] = 360 + 460
K["Tt_3"] = 715 + 460
K["Tt_4"] = 1400 + 460
K["Tt_5"] = 890+ 460
K["Qr"] = 4.326E7
K["cp_3"] = 1.005E3
K["cp_4"] = 1.089E3
K["gamma"] = 1.4
K["R"] = 287
K["mdot_fan"] = 120.202
K["mdot_core"] = 88.451
K["Fn"] = 8.007E4

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
