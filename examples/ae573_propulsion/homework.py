"""AE 573 turbojet cycle homework (``HomeWork``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.

Equations that do not parse (typos in the original) are left out: ``pt_7 =
pi_ AB * pt_5``
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'Dram = mdot_0 * V_0',
    'Fn = (mdot_0 + mdot_f) * V_9 - mdot_0 * V_0 + (p_9-p_0)*A_9',
    'Fg = (mdot_0+mdot_f)*V_9 - mdot_0 * V_0',
    'Fn = Fg - Dram',
    'TSFC = mdot_f / Fn',
    'TSFC = f / (Fn/mdot_0)',
    'eta_th = DeltaKEdot / pth',
    'pth = mdot_f * Qr + mdot_f_AB * Qr_AB',
    'DeltaKEdot = (mdot_0 + mdot_f + mdot_f_AB) * V_9^2 / 2 - mdot_0 * V_0^2 /2',
    'eta_p = Fn * V_0 / DeltaKEdot',
    'eta_p = 2/(1+V_9/V_0)',
    'f_AB = (1+f)*(ht_7-ht_5) / (Qr_AB * eta_AB - ht_7)',
    'eta_p = (Fn/mdot_0)*V0 / ((1+f+f_AB)*V_9^2/2 - v_0^2/2)',
    'Is = Fn/(mdot_f*g_0)',
    'eta_o = (Fn/mdot_0)*V_0 / (f*Qr)',
    'eta_o = (V_0/Qr) / TSFC',
    'pi_n = ((pi_AB*pi_t*pi_b*pi_c*pi_d*pi_r*p_0/p_9)^((gamma-1)/gamma) - eta_n * ((pi_AB*pi_t*pi_b*pi_c*pi_d*pi_r*p_0/p_9)^((gamma-1)/gamma) - 1))^(-gamma/(gamma-1))',
    'pi_d = pt_0 / pt_2',
    'eta_d = ((pt_2/p_0)^((gamma_c-1)/gamma_c) -1) / (M_0^2*(gamma_c-1)/2)',
    'pi_c = pt_3/pt_2',
    'pi_b = pt_4 / pt_3',
    'pi_AB = pt_7 / pt_5',
    'pi_n = pt_9 / pt_7',
    'pi_t = pt_5 / pt_4',
    'pi_r = pt_0 / p_0',
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
    'tau_t = Tt_5 / Tt_4',
    'pi_c_optimum = ((tau_r + tau_lambda) / (2*tau_r))^(gamma/(gamma-1))',
    'eta_d = (ht_2s/h_0 - 1) / (ht_2/h_0 - 1)',
    'eta_d = (Tt_2s/T_0 - 1) / (ht_0/h_0 - 1)',
    'eta_d = ((pt_2/p_0)^((gamma_c-1)/gamma_c) - 1) / (tau_r - 1)',
    'tau_c = pi_c^((gamma-1)/(gamma*e_c))',
    'P_r = mdot_t * (ht_2-ht_3)',
    'F= Fn',
    'gamma_0 = gamma_c',
    'gamma_1 = gamma_c',
    'gamma_2 = gamma_c',
    'gamma_3 = gamma_c',
    'gamma_4 = gamma_t',
    'gamma_5 = gamma_t',
    'gamma_6 = gamma_AB',
    'gamma_7 = gamma_AB',
    'gamma_8 = gamma_AB',
    'gamma_9 = gamma_AB',
    'cp_0 = cp_c',
    'cp_1 = cp_c',
    'cp_2 = cp_c',
    'cp_3 = cp_c',
    'cp_4 = cp_t',
    'cp_5 = cp_t',
    'cp_6 = cp_t',
    'cp_7 = cp_AB',
    'cp_8 = cp_AB',
    'cp_9 = cp_AB',
    'tau_d = 1',
    'tau_n = 1',
    'M_2 = ( ((gamma_c-1)*M_0^2+2) / (2*gamma_c*M_0^2-(gamma_c-1)) )^.5',
    'rho_2 / rho_0 = ((gamma_c+1)*M_0^2) / (2 + (gamma_c-1)*M_0^2)',
    'p_2 / p_0 = 1 + 2*gamma_c/(gamma_c+1) * (M_0^2-1)',
    'T_2 / T_0 = (p_2 / p_0) / (rho_2 / rho_0)',
    'pt_2 / pt_0 = (((gamma_c+1)*M_0^2) / ((gamma_c-1)*M_0^2+2))^(gamma_c/(gamma_c-1)) * ((gamma_c+1) / (2*gamma_c*M_0^2-(gamma_c-1)))^(1/(gamma_c-1))',
    'Tt2 = Tt0',
]

# (templates, token, values): each template is repeated with ``token`` -> ``_<value>``
REPEATED = [
    (
        [
            'M_st = V_st / a_st',
            'a_st = (gamma_st*R*T_st)^.5',
            'pt_st / p_st = (1 + (gamma_st-1)/2 * M_st^2)^(gamma_st/(gamma_st-1))',
            'rhot_st / rho_st = (1 + (gamma_st-1)/2 * M_st^2)^(1/(gamma_st-1))',
            'Tt_st / T_st = (1 + (gamma_st-1)/2 * M_st^2)^1',
            'p_st = rho_st * R * T_st',
        ],
        '_st',
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    ),
]

K = {}
K["mdot_0"] = 100
K["M_0"] = 2
K["p_0"] = 20e3
K["T_0"] = 228
K["gamma"] = 1.4
K["R"] = 287
K["mdot_f"] = 2.8
K["Qr"] = 42000e3
K["V_9"] = 1200
K["p_9"] = K["p_0"]
K["eta_b"] = .995
K["pi_b"] = .95
K["pi_d"] = .9
K["pi_c"] = 12
K["e_c"] = .9
K["g_0"] = 9.81
K["A_9"] = 1
K["gamma_c"] = 1.4
K["gamma_t"] = 1.33
K["gamma_AB"] = 1.3
K["cp_c"] = 1004
K["cp_t"] = 1156
K["cp_AB"] = 1234
K["mdot_f_AB"] = 0
K["Qr_AB"] = 0

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
