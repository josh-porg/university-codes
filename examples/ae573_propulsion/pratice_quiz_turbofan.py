"""AE 573 practice quiz (turbofan) (``Pratice_quiz_turbofan``).

Converted from the MATLAB equation lists: every relation with a single
unknown is solved in turn (``solveRelations``) until nothing changes, then
all solved quantities are printed.
"""

import numpy as np  # noqa: F401

from unicodes.atmosphere import isa
from unicodes.thermo import build_relations, expand_stations, solve_relations

EQUATIONS = [
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
    'cp_c = cv_c + R',
    'cp_t = cv_t + R',
    'cp_c = gamma_t * R / (gamma_t - 1)',
    'cv_c = R / (gamma_t - 1)',
    'cp_t = gamma_t * R / (gamma_t - 1)',
    'cv_t = R / (gamma_t - 1)',
]

# (templates, token, values): each template is repeated with ``token`` -> ``_<value>``
REPEATED = [
    (
        [
            'tau_x = pi_x^ ((gamma_c-1)/(gamma_c *e_x ))',
        ],
        '_x',
        ['d', 'f', 'c', 'fn'],
    ),
    (
        [
            'tau_x = pi_x^ ((gamma_t-1)/(gamma_t *e_x ))',
        ],
        '_x',
        ['t_hp', 't_lp', 'n'],
    ),
    (
        [
            'M_st = V_st / a_st',
            'a_st = (gamma_c*R*T_st)^.5',
            'pt_st / p_st = (1 + (gamma_c-1)/2 * M_st^2)^(gamma_c/(gamma_c-1))',
            'rhot_st / rho_st = (1 + (gamma_c-1)/2 * M_st^2)^(1/(gamma_c-1))',
            'Tt_st / T_st = (1 + (gamma_c-1)/2 * M_st^2)^1',
            'p_st = rho_st * R * T_st',
            'h_st = u_st + p_st / rho_st',
            'h_st = cp_c*Tt_st',
            'p_st = rho_st*a_st^2/gamma_c',
            'ht_st = h_st + V_st^2 / 2',
        ],
        '_st',
        [0, 1, 2, 3, 13, 19],
    ),
    (
        [
            'M_st = V_st / a_st',
            'a_st = (gamma_t*R*T_st)^.5',
            'pt_st / p_st = (1 + (gamma_t-1)/2 * M_st^2)^(gamma_t/(gamma_t-1))',
            'rhot_st / rho_st = (1 + (gamma_t-1)/2 * M_st^2)^(1/(gamma_t-1))',
            'Tt_st / T_st = (1 + (gamma_t-1)/2 * M_st^2)^1',
            'p_st = rho_st * R * T_st',
            'h_st = u_st + p_st /rho_st',
            'h_st = cp_t*Tt_st',
            'p_st = rho_st*a_st^2/gamma_t',
            'ht_st = h_st + V_st^2 / 2',
        ],
        '_st',
        [45, 4, 5, 6, 7, 8, 9],
    ),
]

altitude = 11.2776E3
_atm = isa(altitude)  # getSatndardAtmosphericValues
T, p, rho, c = _atm.T, _atm.p, _atm.rho, _atm.a
K = {}
K["pi_c"] = 40
K["pi_d"] = .995
K["pi_f"] = 1.6
K["pi_fn"] = .98
K["pi_b"] = .95
K["pi_n"] = .98
K["e_c"] = .9
K["e_f"] = .9
K["e_t"] = .9
K["tau_lambda"] = 7
K["eta_b"] = .98
K["eta_m"] = .975
K["M_9"] = 1
K["M_19"] = 1
K["M_0"] = .8
K["T_0"] = T
K["p_0"] = p
K["rho_0"] = rho
K["a_0"] = c
K["gamma_c"] = 1.4
K["gamma_t"] = 1.33
K["cp_c"] = 1.004
K["cp_t"] = 1.156
K["R"] = 287
K["alpha"] = 6
K["Qr"] = 42800

equations = EQUATIONS + [e for t, token, vals in REPEATED for e in expand_stations(t, vals, token)]
relations = build_relations(equations)
solved = solve_relations(relations, K)
for name in sorted(solved, key=str.lower):
    print(f"{name:24s} {solved[name]:.6g}")
