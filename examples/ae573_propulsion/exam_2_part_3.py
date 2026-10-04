"""AE 573 exam 2 part 3 (``Exam_2_part_3``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.

Equations that do not parse (typos in the original) are left out: ``pt_7 =
pi_ AB * pt_5``
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'eta_gb = P_prop/P_lpt*eta_gb',
    'eta_prop=F_prop*V_0/P_prop',
    'eta_m_hpt = (ht_3-ht_2)/((1+f)*(ht_4-ht45))',
    'eta_n = (ht_5-h_5)/(ht5-h5s)',
    'Tt_2 = Tt_0',
    'pt_45/P_t4 = (Tt_45/Tt_4)^(gamma/((gamma-1)*e_t_hpt))',
    'h_8 = h_9',
    'T_8 = T_9',
    'p_8 = p_9',
    'v_8 = v_9',
    'pi_n = pt_9/pt_5',
    'eta_th_ideal = 1-(PR)^(-(gamma-1)/gamma)',
    'P_c = mdot_9*v_9^2/2-mdot_0*v_0^2/2',
    'gamma_c = gamma',
    'gamma_t = gamma',
    'cp_c = cp',
    'cp_t = cp',
    'TSFC = mdot_f/P_prop+P_core',
    'PSFC = mdot_f/(F_prop+F_core)',
    'F_total = F_prop+F_core',
    'tau_c = pi_c^((gamma-1)/(e_c*gamma_c))',
    'Tt_45 = ht_45/cp',
    'eta_m_hpt*(1+f)*(ht_4-ht_45) = ht_3-ht_2',
    'pt_4 = pt_3*eta_b',
    'pi_n = pt_9/pt_5',
    'v_9^2 = 2 * cp *(Tt_5-T_9)',
    'pi_t = tau_t^(gamma/e_t*(gamma-1))',
    'P_fan = P_prop',
    'P/(mdot_0*(1+f)) = ht_45*(1-(p_0/pt_45)^((gamma-1)/gamma))',
    'P/(mdot_0*(1+f)) = P_i_tot/mdot_9',
    'ht_45*(1-(p_0/pt_45)^((gamma-1)/gamma)) = P_i_tot/mdot_9',
    'mdot_0 = mdot_9/(1+f)',
    'f = (ht_4-ht_3)/(Qr*eta_b-ht_4)',
    'alpha = ((P_lpt/mdot_9)/eta_lpt)/(P_i_tot/mdot_9)',
    'P_lpt = mdot9*(ht_45-ht_5)',
    'P_lpt = mdot_9*eta_lpt*alpha*ht_45*(1-(p_0/pt_45)^(gamma-1)/gamma)',
    'P_lpt = P_prop*eta_prop*eta_gb',
    'P_prop = eta_gb*eta_m_lpt*P_lpt',
    'P_prop = mdot_0 *(1+f)*eta_gb*eta_m_lpt*eta_lpt*alpha*ht_45*(1-(p_0/pt_45)^(gamma-1)/gamma)',
    'P_core = mdot_0/2*((1+f)*v_9^2-v_0^2)',
    'ht_45 = tau_hpt*tau_lambda*h_0',
    'F_prop*V_0 = eta_prop*eta_gb*eta_m_lpt*P_lpt',
    'F_core = mdot_0*((1+f)*v_9-v_0)',
    'P_prop/eta_gb = F_prop*v_0',
    'v_9 = 2*(1-alpha)*eta_n*ht_45*(1-(p_9/pt_45)^((gamma-1)/gamma))',
    'mdot_9 = mdot_0 * (1+f)',
    'p_0 = p_9',
    'eta_o = eta_th*eta_p',
    'eta_p = F_total*v_0/(P_prop+P_core)',
    'eta_th = (P_prop+P_core)/(mdot_f*Qr)',
    'f = (ht_4-ht_3)/(Qr*eta_b-ht_4)',
    'NPR = pt_7/p_0',
    'pi_n = pt_9/pt_7',
    'NPR_crit = 1/pi_cn*((gamma+1)/2)^(gamma/(gamma-1))',
    'NPR_crit = pt_7/p_0',
    'M_9 = (2/(gamma-1)*(pt_9/p_9)^((gamma-1)/gamma)-1)',
    'p8rat = (1 + (gamma-1)/2 * M_8^2)^(gamma/(gamma-1))*((gamma+1)/2)^(gamma/(gamma-1))',
    'p8rat = p8/p0',
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
    'pi_b = pt_4/pt_3',
    'tau_c = Tt_3 / Tt_2',
    'tau_f = Tt_2 / Tt_1',
    'tau_f = Tt_13 / Tt_1',
    'cp = cv + R',
    'cp = gamma * R / (gamma - 1)',
    'cv = R / (gamma - 1)',
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
]

# (templates, token, values): each template is repeated with ``token`` -> ``_<value>``
REPEATED = [
    (
        [
            'tau_x = pi_x^ ((gamma -1)/(gamma *e_x ))',
            'v_x = V_x',
        ],
        '_x',
        ['d', 'f', 'c', 'fn', 'n'],
    ),
    (
        [
            'v_st = V_st',
            'M_st = V_st / a_st',
            'a_st = (gamma*R*T_st)^.5',
            'pt_st / p_st = (1 + (gamma-1)/2 * M_st^2)^(gamma/(gamma-1))',
            'rhot_st / rho_st = (1 + (gamma-1)/2 * M_st^2)^(1/(gamma-1))',
            'Tt_st / T_st = (1 + (gamma-1)/2 * M_st^2)^1',
            'ht_st = Tt_st * cp',
            'p_st = rho_st * R * T_st',
            'V_9 = a_0*M_9*(tau_lambda*tau_t/(1+(gamma_c-1)/2*M_9^2))^.5',
            'V_19 = a_0*M_19*(tau_lambda*tau_t/(1+(gamma_c-1)/2*M_19^2))^.5',
            'M_9 = (2/(gamma_t-1)*(pi_n*pi_t*pi_b*pi_c*pi_d*pi_r*p_0/p_9)^(gamma_t/(gamma_t-1))-1)^.5',
            'M_19 = (2/(gamma_t-1)*(pi_fn*pi_f*pi_b*pi_r*p_0/p_19)^(gamma_c/(gamma_c-1))-1)^.5',
            'P_c/P_f = (tau_c - 1) / (alpha*(tau_f-1))',
            'f = (tau_lambda - tau_c*tau_r) / (Qr*eta_b/(cp_t*T_0)-tau_lambda)',
            'eta_m * mdot_0 * (1+f) * (h_t4-ht5) = mdot_0 * (ht_3-ht_2) + alpha * mdot_0 * (ht_13-ht_2)',
            'eta_m * ht_4/h_0 * (1+f) * (1-tau_t) = ht_2 / h_0 * ( (tau_c-1) + alpha * (tau_f-1))',
            'tau_t = 1 - tau_r * ((tau_r-1) + alpha*(tau_f-1)) / (eta_m*(1+f)*tau_lambda)',
            'tau_lambda = ht_4/h_0',
            '(tau_lambda)^.5 = tau_r * tau_c',
            'tau_r = Tt_0/T_0',
            'pi_r = Tt_2/Tt_0',
            'tau_t = tau_t_hp*tau_t_lp',
            'f = (tau_lambda - tau_r*tau_c) / (Qr/h_0 * eta_b - tau_lambda)',
            'Tt_1 = Tt_0',
            'p_19 = p_0',
            'alpha = mdot_f / mdot_c',
            'tau_c = pi_c^ ((gamma_t-1)/(gamma_t*e_c))',
            'tau_f = pi_f^ ((gamma_t-1)/(gamma_t*e_f))',
            'P_f = mdot_f *cp_c * (Tt_1 - Tt_13)',
            'P_c = mdot_c *cp_c * (Tt_2 - Tt_3)',
            'mdot_c = mdot_2',
            'mdot_f = mdot_2',
            'pi_d = pt_1/pt_0',
            'pi_f = pt_13/pt_1',
            'pi_fn = pt_19/pt_13',
            'pi_c = pt_3/pt_2',
            'tau_c = Tt_3 / Tt_2',
            'tau_f = Tt_2 / Tt_1',
            'tau_f = Tt_13 / Tt_1',
        ],
        '_st',
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 19],
    ),
]

K = {}
K["pi_c"] = 25
K["pi_d"] = .98
K["pi_f"] = 1.5
K["eta_gb"] = .996
K["eta_prop"] = .75
K["e_c"] = .9
K["M_0"] = .76
K["p_0"] = 35E3
K["T_0"] = 242
K["gamma"] = 1.4
K["cp"] = 1004
K["R"] = 287
K["Tt_4"] = 1600
K["Qr"] = 42000
K["eta_b"] = .99
K["pi_b"] = .95
K["eta_m_hpt"] = .99
K["e_hpt"] = .85
K["alpha"] = .87
K["eta_lpt"] = .9
K["eta_m_lpt"] = .985
K["eta_n"] = .96
K["pi_t"] = 1
K["pi_n"] = 1

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
