"""AE 573 quiz 2 (``Quiz_2``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'pi_n = ((pi_AB*pi_t*pi_b*pi_c*pi_d*pi_r*p_0/p_9)^((gamma-1)/gamma) - eta_n * ((pi_AB*pi_t*pi_b*pi_c*pi_d*pi_r*p_0/p_9)^((gamma-1)/gamma) - 1))^(-gamma/(gamma-1))',
    'pi_d = pt_0 / pt_2',
    'eta_d = ((pt_2/p_0)^((gamma_c-1)/gamma_c) -1) / (M_0^2*(gamma_c-1)/2)',
    'pi_c = pt_3/pt_2',
    'pi_c_lp = pt_2_5/+pt_2',
    'pi_c_hp = pt_3/pt_2_5',
    'pi_c = pi_c_lp*pi_c_hp',
    'pi_b = pt_4 / pt_3',
    'pi_AB = pt_7 / pt_5',
    'pi_n = pt_9 / pt_7',
    'pi_t = pt_5 / pt_4',
    'pi_r = pt_0 / P_0',
    'pi_r = (1+(gamma_c-1)/2*M_0^2)^(gamma_c/(gamma_c-1))',
    'tau_r = Tt_0 / T_0',
    'tau_r = 1 + (gamma_c-1)/2 *M_0^2',
    'tau_d = Tt_2/Tt_0',
    'tau_lambda = cp_4*Tt_4 / cp_0*T_0',
    'eta_b = Qr_actual / Qr_ideal',
    'eta_AB = Qr_AB_actual / Qr_AB_ideal',
    'tau_lambda_AB = cp_AB*Tt_7 / cp_0*T_0',
    'tau_n = Tt_9 / Tt_7',
    'tau_c = Tt_3 / Tt_2',
    'tau_c = tau_c_lp*tau_c_hp',
    'tau_c_lp = Tt_2 / Tt_2_5',
    'tau_c_hp = Tt_3 / Tt_2_5',
    'tau_t = Tt_5 / Tt_4',
    'TSFC = mdot_f / Fn',
    'TSFC = f / (Fn/mdot_0)',
    'eta_th = DeltaKEdot / pth',
    'pth = mdot_f * Qr',
    'DeltaKEdot = (mdot_0 + mdot_f) * V_9^2 / 2 - mdot_0 * V_0^2 /2',
    'eta_p = Fn * V_0 / DeltaKEdot',
    'eta_p = 2/(1+V_9/V_0)',
    '(1+f)*cp_4*Tt_4 - cp_3*Tt_3 = f * Qr',
    'Tt_1 = Tt_0',
    'Dram = mdot_0 * V_0',
    'Fn = (mdot_0 + mdot_f) * V_9 - mdot_0 * V_0 + (p_9-p_0)*A_9',
    'Fg = (mdot_0+mdot_f)*V_9 - mdot_0 * V_0',
    'Fn = Fg - Dram',
    'TSFC = mdot_f / Fn',
    'Tt_1 = Tt_0',
    'p_9 = p_0',
    'tau_t = pi_t^ ((gamma_t-1)/(gamma_t)*e_t)',
    'tau_c = pi_c^ ((gamma_c-1)/(gamma_c*e_c))',
    'P_compressor = mdot_core *cp * (Tt_2 - Tt_3)',
    'mdot_core = mdot_2',
    'pi_d = pt_1/pt_0',
    'pi_c = pt_3/pt_2',
    'tau_c = Tt_3 / Tt_2',
    'tau_5 = Tt_4 / Tt_5',
    'gamma_0 = gamma_c',
    'gamma_1 = gamma_c',
    'gamma_2 = gamma_c',
    'gamma_3 = gamma_c',
    'gamma_4 = gamma_t',
    'gamma_5 = gamma_t',
    'gamma_6 = gamma_t',
    'gamma_7 = gamma_t',
    'gamma_8 = gamma_t',
    'gamma_9 = gamma_t',
    'gamma_d = gamma_c',
    'gamma_b = gamma_t',
    'gamma_t = gamma_t',
    'gamma_n = gamma_t',
    'gamma_AB = gamma_t',
    'cp_0 = cp_c',
    'cp_1 = cp_c',
    'cp_2 = cp_c',
    'cp_3 = cp_c',
    'cp_4 = cp_t',
    'cp_5 = cp_t',
    'cp_6 = cp_t',
    'cp_7 = cp_t',
    'cp_8 = cp_t',
    'cp_9 = cp_t',
    'cp_d = cp_c',
    'cp_b = cp_t',
    'cp_t = cp_t',
    'cp_n = cp_t',
    'cp_AB = cp_t',
]

# (templates, token, values): each template is repeated with ``token`` -> ``_<value>``
REPEATED = [
    (
        [
            'cp_x = cv_x + R',
            'cp_x = gamma_x * R / (gamma_x - 1)',
            'cv_x = R / (gamma_x - 1)',
        ],
        '_x',
        ['d', 'c', 'b', 't', 'AB', 'n'],
    ),
    (
        [
            'M_st = V_st / a_st',
            'a_st = (gamma_st*R*T_st)^.5',
            'pt_st / p_st = (1 + (gamma_st-1)/2 * M_st^2)^(gamma_st/(gamma_st-1))',
            'rhot_st / rho_st = (1 + (gamma_st-1)/2 * M_st^2)^(1/(gamma_st-1))',
            'Tt_st / T_st = (1 + (gamma_st-1)/2 * M_st^2)^1',
            'p_st = rho_st * R * T_st',
            'cp_st = cv_st + R',
            'cp_st = gamma_st * R / (gamma_st - 1)',
            'cv_st = R / (gamma_st - 1)',
        ],
        '_st',
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    ),
]

K = {}
K["pi_c"] = 15
K["pi_d"] = .9
K["pi_b"] = .95
K["pi_t"] = .936
K["pi_n"] = .95
K["pi_AB"] = .98
K["e_c"] = .9
K["e_t"] = .85
K["eta_b"] = .99
K["Qr"] = 42000E3
K["M_0"] = 2
K["p_0"] = 25E3
K["T_0"] = 228.15
K["gamma_c"] = 1.4
K["gamma_t"] = 1.33
K["R"] = 287
K["mdot_0"] = 100
K["mdot_core"] = 100
K["M_9"] = 2.62
K["T_9"] = 603
K["Tt_4"] = 1760
K["Fn"] = 4588.86

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
