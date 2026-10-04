"""AE 573 final exam problem 1 (``FinalExam_1``).

mixed after burning turbofan

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.

Equations that do not parse (typos in the original) are left out:
``p_6*(1+gamma_6*M_6^2)*(1+A15/A5) -
p5((1+gamma_t*M_5^2)+(A15/A5)*(1+gamma_c*M_15^2)) = 0``
"""

import numpy as np  # noqa: F401

from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
    'Tt_13 = Tt_15',
    'Dram = mdot_0 * V_0',
    'Fn = (mdot_0 + mdot_f) * V_9 - mdot_0 * V_0 + (p_9-p_0)*A_9',
    'Fg = (mdot_0+mdot_f)*V_9 - mdot_0 * V_0',
    'Fn = Fg - Dram',
    'TSFC = mdot_f / Fn',
    'TSFC = f / (Fn/mdot_0)',
    'eta_th = DeltaKEdot / pth',
    'pth = mdot_f * Qr + mdot_f_AB * Qr_AB',
    'DeltaKEdot = (mdot_0 + mdot_f + mdot_f_AB) * V_9^2 / 2 - mdot_0 * V_0^2 /2',
    'Tt_1 = Tt_0',
    'pt_1 = pt_0',
    'rhot_1 = rhot_0',
    'T_1 = T_0',
    'p_1 = p_0',
    'rho_1 = rho_0',
    'M_0 = M_1',
    'ht_0 = ht_1',
    'h_0 = h_1',
    'V_0 = V_1',
    'a_0 = a_1',
    'p_9 = p_0',
    'alpha = mdot_fan / mdot_core',
    'tau_c = pi_c^ ((gamma_c-1)/(gamma_c*e_c))',
    'tau_f = pi_f^ ((gamma_c-1)/(gamma_c*e_f))',
    'P_fan = mdot_fan *cp * (Tt_1 - Tt_13)',
    'P_compressor = mdot_core *cp * (Tt_2 - Tt_3)',
    '(1+alpha)*mdot_0 = mdot_in_total',
    'mdot_core = mdot_0',
    'mdot_2 = mdot_core + mdot_fan',
    'mdot_1 = mdot_in_total',
    'mdot_6 = mdot_in_total+mdot_f',
    'pi_d = pt_1/pt_0',
    'pi_f = pt_13/pt_1',
    'pi_f = pt_2/pt_1',
    'pi_fd = pt_15/pt_13',
    'pi_c = pt_3/pt_2',
    'pi_inlet_outlet = p_9/p_0',
    'tau_c = Tt_3 / Tt_2',
    'tau_f = Tt_2 / Tt_1',
    'tau_f = Tt_13 / Tt_1',
    'tau_fd = Tt_15 / Tt_13',
    'pt_5 = pt_13',
    'eta_f = (pi_f^((gamma_c-1)/gamma_c)-1)/(tau_f-1)',
    'tau_f = 1 + 1/(eta_f)*(pi_f^((gamma_c-1)/gamma_c)-1)',
    'pi_f = (1+eta_f*(tau_f-1))^(gamma_c/(gamma_c-1))',
    'eta_m*mdot_1*(f+1)*(ht_4-ht_5) = mdot_0*(ht_3-ht_2) + alpha*mdot_0*(ht_13-ht_2)',
    'mdot_15/mdot_5 = alpha/(1+f)',
    'mdot_15 = alpha*mdot_0',
    'mdot_5 = mdot_0+mdot_f',
    'mdot_6 = (1 + alpha)*mdot_0 + mdot_f',
    'p_5 = p_15',
    'mdot_5 *ht_5 + mdot_15*ht_15 = mdot_6 * ht_6',
    'h_t6 = ((1+f)*h_t5 + alpha*ht_15) / (1+alpha+f)',
    'ht_6 = h_0 + ((1+f)*tau_t*tau_lambda + alpha*tau_r*tau_f) / (1+alpha+f)',
    'F= Fn',
    'gamma_0 = gamma_c',
    'gamma_1 = gamma_c',
    'gamma_2 = gamma_c',
    'gamma_3 = gamma_c',
    'gamma_13 = gamma_c',
    'gamma_15 = gamma_c',
    'gamma_4 = gamma_t',
    'gamma_5 = gamma_t',
    'gamma_6 = (mdot_5*cp_t + mdot_15*cp_c)/(mdot_5*(cp_t/gamma_t)+mdot_15*(cp_c/gamma_c))',
    'gamma_7 = gamma_AB',
    'gamma_8 = gamma_AB',
    'gamma_9 = gamma_AB',
    'cp_0 = cp_c',
    'cp_1 = cp_c',
    'cp_2 = cp_c',
    'cp_3 = cp_c',
    'cp_13 = cp_c',
    'cp_15 = cp_c',
    'cp_4 = cp_t',
    'cp_5 = cp_t',
    'cp_6 = (mdot_5*cp_t + mdot_15*cp_c)/(mdot_5+mdot_15)',
    'cp_7 = cp_AB',
    'cp_8 = cp_AB',
    'cp_9 = cp_AB',
    'tau_d = 1',
    'tau_n = 1',
    'Tt_max = Tt_9',
    'Tt_min = Tt_2',
    'pt_0 = p_0 * (1+(gamma_c-1)/2*M_0^2)^(gamma_c/(gamma_c-1))',
    'pi_r = (1+(gamma_c-1)/2*M_0^2)^(gamma_c/(gamma_c-1))',
    'Tt_0 = T_0 * (1+(gamma_c-1)/2*M_0^2)',
    'tau_r = (1+(gamma_c-1)/2*M_0^2)',
    'a_0 = sqrt(gamma_c*R*T_0)',
    'pt_2 = pi_d*pt_0',
    'pt_2 =  p_0 * (1+eta_d*(gamma_c-1)/2*M_0^2)^(gamma_c/(gamma_c-1))',
    'Tt_2 = Tt_0',
    'pt_3 = pt_2 * pi_c',
    'Tt_3 = Tt_2 * pi_c^((gamma_c-1)/(gamma_c*e_c))',
    'tau_c = pi_c^((gamma_c-1)/(gamma_c*e_c))',
    'pt_4 = pt_3 * pi_b',
    'ht_4 = cp_t * Tt_4',
    'ht_4 = 1/(f+1) * (ht_3 + f*Qr*eta_b)',
    'f = (ht_4-ht_3)/(Qr*eta_b-ht_4)',
    'f = (tau_lambda-tau_r*tau_c)/(Qr*eta_b/h_0-tau_lambda)',
    'tau_lambda = ht_4/h_0',
    'pt_13 = pt_2 * pi_f',
    'Tt_13 = Tt_2 * pi_f^((gamma_c-1)/(gamma_c*e_f))',
    'tau_f = pi_f^((gamma_c-1)/(gamma_c*e_f))',
    'pt_15 = pt_13 * pi_fd',
    'Tt_15 = Tt_13',
    'ht_15 = ht_13',
    'tau_fd = 1',
    'et_m*(f+1)*(ht_4-ht_5)  = (ht_3-ht_2) + alpha*(ht_13-ht_2)',
    'eta_m*(1+f)*tau_lambda*(1-tau_t) = tau_r*(tau_c-1) + alpha*tau_r*(tau_f-1)',
    'pt_5/pt_2 = pi_t * pi_b * pi_c',
    'pt_15/pt_2 = pi_df * pi_f',
    'pi_t = pi_fd*pi_f/(pi_b*pi_c)',
    'tau_t = pi_t^((gamma_t-1)*e_t/(gamma_t))',
    'tau_t = (pi_fd*pi_f/(pi_b*pi_c))^((gamma_t-1)*e_t/(gamma_t))',
    'alpha = (eta_m*(1+f)*tau_lambda*(1-tau_t)-tau_r*(tau_c-1)) / tau_r*(tau_f-1)',
    'ht_6 = h_0 * ((1+f)*tau_t*tau_lambda+alpha*tau_f*tau_r) / (1+alpha+f)',
    'pt_15 = p_15 * (1+(gamma_c-1)/2*M_15^2)^(gamma_c/(gamma_c-1))',
    'pt_15 = p_5 * (1+(gamma_c-1)/2*M_5^2)^(gamma_c/(gamma_c-1))',
    'M_15^2 = 1/(gamma_c-1) * ((1+(gamma_t-1)/2*M_5^2)^((gamma_c-1)/gamma_c)-1)',
    'mdot_15 = gamma_c *p_15*A15*M_15/a_15',
    'alpha*mdot_0 = gamma_c *p_15*A15*M_15/a_15',
    'mdot_15 = alpha*mdot_0',
    'mdot_5 = gamma_t *p_5*A5*M_5/a_5',
    '(1+f) * mdot_0 = gamma_t *p_5*A5*M_5/a_5',
    'mdot_5 = (1+f)*mdot_0',
    'A_15/A_5 = alpha/(1+f)*(gamma_t/gamma_c)*(a_15/a_5)*(M_5/M_15)',
    '(gamma_6*p_6*M_6*(1+(gamma_6-1)*M_6^2/2)^(1/2)) / (((gamma_6-1)*cp_6*Tt_6)^(1/2)) * (1+A15/A5) = gamma_t*p_5*M_5/a_5 + gamma_c*p_15*M_15/a_15',
    'p_6/p_5 * (1+gamma_6*M_6^2) = C1',
    '((1+gamma_t*M_5^2)+(A15/A5)*(1+gamma_c*M_15^2)) / (1+A15/A5) = C1',
    'p_6/p_5 * M_6*(1+(gamma_6-1)*M_6^2/2)^(1/2) = C2',
    '((gamma_t/gamma_6)*(M_5/a_5)+(gamma_c/gamma_6)*(M_15*(A15/A5)/a_15)) * ((gamma_6-1)*cp_6*Tt_6)^(1/2) = C2',
    '(1+gamma_6*M_6^2) / (M_6*(1+(gamma_6-1)*M_6^2/2)^(1/2)) = (C)^(1/2)',
    'C1/C2 = (C)^(1/2)',
    '(1+gamma_6*M_6^2)^2 = C * M_6^2 * (1+(gamma_6-1)*M_6^2/2)',
    '(gamma_6^2-C*(gamma_6-1)/2)*M_6^4 + (2*gamma_6-C)*M_6^2+1 = 0',
    'M_6^2 = (C - 2*gamma_6 - ((C-2*gamma_6)^2 - 4 *(gamma_6^2-C*(gamma_6-1)/2))^(1/2)) / (2*gamma_6^2-C*(gamma_6-1))',
    'p_6/p_5 = C1/(1+gamma_6*M_6^2)',
    'pi_m_ideal = pt_6/pt_5',
    'pi_m_ideal = p_6/p_5 * (1+(gamma_6-1)*M_6^2/2)^(gamma_6/(gamma_6-1)) / (1+(gamma_t-1)*M_5/2)^(gamma_t/(gamma_t-1))',
    'pi_m = pi_m_ideal * pi_m_frictionalLoss',
    '(mdot_6 + mdot_f_AB)*ht_7 - mdot_6*ht_6 = mdot_f_AB * Qr_AB * eta_AB',
    'f_AB = mdot_f_AB / mdot_6',
    'f_AB = (ht_7 - ht_6) / (Qr_AB*eta_AB-ht_7)',
    'f_AB = (tau_lambda_AB - ht_6/h_0) / (Qr_AB*eta_AB/h_0 - tau_lambda_AB)',
    'pt_7 = pt_6*pi_AB',
    'Tt_9 = Tt_7',
    'pt_9 = pt_7 * pi_n',
    'eta_n = ((pt_7/p_9)^((gamma_AB-1)/gamma_AB) - pi_n^((gamma_AB-1)/-gamma_AB)) / ((pt_7/p_9)^((gamma_AB-1)/gamma_AB) - 1)',
    'M_9^2 = 2/(gamma_AB - 1) * ((pt_9/p_9)^((gamma_AB-1)/gamma_AB) - 1)',
    'a_9^2 = gamma_AB*R_AB*Tt_9 / (1 + (gamma_AB-1)*M_9/2)',
    'V_9 = M_9*a_9',
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
            'ht_st = h_st + V_st^2/2',
            'st_st = s_st',
            'Tt_st = T_st + V_st^2/(2*cp_st)',
            'p_st = rho_st * R * T_st',
            'cp_st = cv_st + R',
            'cp_st = gamma_st * R / (gamma_st - 1)',
            'cv_st = R / (gamma_st - 1)',
            'ht_st = cp_st*ht_st',
        ],
        '_st',
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 15],
    ),
]

K = {}
K["f"] = .013
K["mdot_6"] = 25+25*.013
K["Tt_5"] = 1.673e3
K["V_9"] = 1.09763e3
K["M_0"] = 2.2
K["p_0"] = 10e3
K["T_0"] = 223.15
K["R"] = 287
K["gamma"] = 1.4
K["gamma_c"] = 1.4
K["gamma_t"] = 1.4
K["gamma_AB"] = 1.4
K["pi_d"] = .85
K["mdot_in_total"] = 25
K["pi_c"] = 15
K["e_c"] = .9
K["pi_f"] = 1.9
K["e_f"] = .9
K["pi_b"] = .95
K["eta_b"] = .98
K["Qr"] = 42800e3
K["Tt_4"] = 1.673e3
K["e_t"] = .92
K["eta_m"] = .95
K["M_5"] = .5
K["pi_m_frictionalLoss"] = .98
K["Qr_AB"] = 0
K["pi_AB"] = 1
K["eta_AB"] = 1
K["mdot_f_AB"] = 0
K["R_AB"] = K["R"]
K["pi_n"] = 1
K["pi_inlet_outlet"] = 1
K["pi_m"] = 1

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
