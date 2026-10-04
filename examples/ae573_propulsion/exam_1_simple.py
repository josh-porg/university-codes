"""AE 573 exam 1 (simplified equation set) (``EXAM_1_simple``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'f = (tau_lambda - tau_c*tau_r) / (Qr*eta_b/(cp*T_0)-tau_lambda)',
    'M_9 = (2/(gamma-1)*(pi_n*pi_t*pi_b*pi_c*pi_d*pi_r*p_0/p_9)^(gamma/(gamma-1))-1)^.5',
    '1-tau_t = P_turbine / (mdot_core * cp * Tt_4)',
    'P_fan + P_compressor = eta_m * P_turbine',
    'Fg_fan = (mdot_fan)*V_19 - mdot_fan * V_0',
    'Tt_1 = Tt_0',
    'p_9 = p_0',
    'p_19 = p_0',
    'alpha = mdot_fan / mdot_core',
    'tau_c = pi_c^ ((gamma-1)/(gamma*e_c))',
    'tau_f = pi_f^ ((gamma-1)/(gamma*e_f))',
    'tau_t = pi_t^ ((gamma-1)/(gamma*e_t))',
    'P_fan = mdot_fan *cp * (Tt_1 - Tt_13)',
    'P_compressor = mdot_core *cp * (Tt_2 - Tt_3)',
    'mdot_core = mdot_2',
    'mdot_fan = mdot_2',
    'pi_d = pt_1/pt_0',
    'pi_f = pt_13/pt_1',
    'pi_f = pt_2/pt_1',
    'pi_fn = pt_19/pt_13',
    'pi_c = pt_3/pt_2',
    'pi_b = pt_4/pt_3',
    'pi_t = pt_5/pt_4',
    'pi_n = pt_9/pt_5',
    'pi_r = pt_0 / p_0',
    'tau_c = Tt_3 / Tt_2',
    'tau_f = Tt_2 / Tt_1',
    'tau_f = Tt_13 / Tt_1',
    'tau_fn = Tt_19 / Tt_13',
    'tau_b = Tt_4/Tt_3',
    'tau_t = Tt_5/Tt_4',
    'tau_n = Tt_9/Tt_5',
    'tau_r = Tt_0 / T_0',
    'tau_lambda = Tt_4/Tt_0',
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
        ['d', 'f', 'c', 'fn', 'n', 't', 'b'],
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
        [0, 1, 2, 3, 4, 5, 13, 19, 9],
    ),
]

K = {}
K["pi_c"] = 45
K["pi_d"] = .99
K["pi_f"] = 1.4
K["pi_fn"] = .98
K["pi_b"] = .97
K["pi_n"] = .98
K["e_c"] = .9
K["e_f"] = .9
K["e_t"] = 1/.8
K["eta_b"] = .99
K["eta_m"] = .99
K["M_0"] = .86
K["p_0"] = 30E3
K["T_0"] = 248
K["gamma"] = 1.4
K["gamma"] = 1.33
K["R"] = 287
K["M_19"] = 1
K["Tt_4"] = 1625
K["Qr_4"] = 42000
K["mdot_fan"] = 240
K["mdot_core"] = 30
K["Qr"] = 42000E3
K["e_fn"] = 1

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
