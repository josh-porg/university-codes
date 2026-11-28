clear
clc
%{
AE 571 Final Project
Group Members: Maggie Bonham, Evelyn Horst, Abha Kaushik, Ben Svoboda,
Jacob Stockley, Jeremy Wegiel

Variable Name Key:
e = Compression Ratio
ER = Equivalence Ratio
Tn = Temperature at state n e.g. T1 is the temperature at state 1
pn = Pressure at state n e.g. P1 is the pressure at state 1
vn = Specific Volume at state n e.g. v1 is the specific volume at state 1
n1 = Polytropic index for compression
n2 = Polytropic index for expansion
alpha = Ratio of pressure increase
beta = Ratio of volume increases
W_net = Net work of the cycle
mep = Mean effective pressure of the cycle
n_th = Thermal efficiency of the cycle
x = Number of carbon atoms in the fuel
y = Number of hydrogen atoms in the fuel
a_stoich = Stoichiometric value of a
a = Number of moles of air in reactants
...

To differentiate between values pertaining to part I and values pertaining
to part II, add a "_1" or "_2" to the end of the variable name e.g. W_net_1
is the net work of the cycle for part I (when gasoline is used).  For
quantities which are constant between parts I and II, such as molecular
weight, nothing is added.

All values are in base SI units unless otherwise specified.
%}

%% Part I
%============================= Question 1 =================================
% Define given values:
e = 16;
ER = 0.8;
T1 = 298;
p1 = 100000;
n1 = 1.38;
n2 = 1.25;
x_1 = 8.26;
y_1 = 15.5;

% Find all values for the chemical equation: 
% CxHy + a(O2 + 3.76N) --> bCO2 + cH2O + dO2 + 3.76 aN2
a_stoich_1 = x_1 + y_1/4;
a_1 = a_stoich_1 / ER;
b_1 = x_1;
c_1 = y_1/2;
d_1 = a_1 - x_1 -(y_1/4);
fprintf('Part I Results:\n')
fprintf('C8H18 + %0.2f*(O2+ 3.76N) --> %0.2f*CO2 + %0.2f*H2O + %0.2f*O2 + %0.2f*N2\n', a_1, b_1, c_1, d_1, 3.76*a_1);


%============================= Question 2 =================================
% First calculate the value of T2_1:
R_air = 287;
R_u = 8.3144598;
v1 = R_air*T1/p1;
v2 = v1/e;
p2 = p1 * e^n1;
T2 = T1 * e^(n1 - 1);

%{
The process for finding the lower heating value is very involved:
hLVH_1 = hbar_1 / MW_fuel_1
hbar_1 = Hc_1 / n_fuel_1
Hc_1 = (HR_1 - HP_1)@T_ref
HR_1 = sum(hbar_species)@T_ref + integral(cvbar_species*dT)
hbar_fuel = 4814*(a1*theta_1 + (1/2)*a2*theta_1^2 + (1/3)*a3*theta_1^3 +
(1/4)*a4*theta_1^4 - a5*theta_1^-1 + a6
hbar_species = Ru*(a6+a1*T + (1/2)*a2*T^2 + (1/3)*a3*T^3 + (1/4)*a4*T^4 +
(1/5)*a5*T^5)
cvbar = cpbar - Ru
cpbar = 4.184*(a1 + (1/2)*a2*theta_1^2 + (1/3)*a3*theta_1^3 +
(1/4)*a4*theta_1^4 - a5*theta_1^-1 + a6)
theta_1 = T2_1/1000
a1 through a6 are found in tables
%}
% Find enthalpy of each species at T2:
theta_1 = T1/1000;
theta_2 = T2/1000;
MW_Gas = 114.8;
a1_Gas = -24.078;
a2_Gas = 256.63;
a3_Gas = -201.68;
a4_Gas = 64.750;
a5_Gas = 0.5808;
a6_Gas = -27.562;
hbar_Gas_T1 = 4184*(a1_Gas + (1/2)*a2_Gas*theta_1^2 + (1/3)*a3_Gas*theta_1^3 + (1/4)*a4_Gas*theta_1^4 - a5_Gas*theta_1^-1 + a6_Gas);

% hbar_Gas_T2 = hbar_Gas_T1 + integral(cpdT)
cpbar_Gas = 4.184*((a1_Gas*theta_2 + (1/2)*a2_Gas*theta_2^2 + (1/3)*a3_Gas*theta_2^3 + (1/4)*a4_Gas*theta_2^4 - a5_Gas*theta_2^-1) - (a1_Gas*theta_1 + (1/2)*a2_Gas*theta_1^2 + (1/3)*a3_Gas*theta_1^3 + (1/4)*a4_Gas*theta_1^4 - a5_Gas*theta_1^-1));
hbar_Gas_T2 = hbar_Gas_T1 + cpbar_Gas;

MW_O2 = 32;
a1_O2_T2 = 0.03212936e+02;
a2_O2_T2 = 0.11274864e-02;
a3_O2_T2 = -0.05756150e-05;
a4_O2_T2 = 0.13138773e-08;
a5_O2_T2 = -0.08768554e-11;
a6_O2_T2 = -0.10052490e04;
hbar_O2_T2 = R_u*(a6_O2_T2 + a1_O2_T2*T2 + (1/2)*a2_O2_T2*T2^2 + (1/3)*a3_O2_T2*T2^3 + (1/4)*a4_O2_T2*T2^4 + (1/5)*a5_O2_T2*T2^5);

MW_N2 = 28;
a1_N2_T2 = 0.03298677e02;
a2_N2_T2 = 0.140822404e-02;
a3_N2_T2 = -0.03963222e-04;
a4_N2_T2 = 0.05641515e-07;
a5_N2_T2 = -0.02444854e-10;
a6_N2_T2 = -0.10208999e04;
hbar_N2_T2 = R_u*(a6_N2_T2 + a1_N2_T2*T2 + (1/2)*a2_N2_T2*T2^2 + (1/3)*a3_N2_T2*T2^3 + (1/4)*a4_N2_T2*T2^4 + (1/5)*a5_N2_T2*T2^5);

MW_CO2 = 44;
a1_CO2_T2 = 0.02275724e02;
a2_CO2_T2 = 0.09922072e-01;
a3_CO2_T2 = -0.10409113e-04;
a4_CO2_T2 = 0.06866686e-07;
a5_CO2_T2 = -0.02117280e-10;
a6_CO2_T2 = -0.04837314e06;
hbar_CO2_T2 = R_u*(a6_CO2_T2 + a1_CO2_T2*T2 + (1/2)*a2_CO2_T2*T2^2 + (1/3)*a3_CO2_T2*T2^3 + (1/4)*a4_CO2_T2*T2^4 + (1/5)*a5_CO2_T2*T2^5);

MW_H2O = 18;
a1_H2O_T2 = 0.03386842e02;
a2_H2O_T2 = 0.03474982e-01;
a3_H2O_T2 = -0.06354696e-04;
a4_H2O_T2 = 0.06968581e-07;
a5_H2O_T2 = -0.02506588e-10;
a6_H2O_T2 = -0.03020811e06;
hbar_H2O_T2 = R_u*(a6_H2O_T2 + a1_H2O_T2*T2 + (1/2)*a2_H2O_T2*T2^2 + (1/3)*a3_H2O_T2*T2^3 + (1/4)*a4_H2O_T2*T2^4 + (1/5)*a5_H2O_T2*T2^5);

% Now combine all the individual enthalpies:
HR_1 = hbar_Gas_T2 + a_1*hbar_O2_T2 + a_1*3.76*hbar_N2_T2;
HP_1 = x_1*hbar_CO2_T2 + (1/2)*y_1*hbar_H2O_T2 + a_1*3.76*hbar_N2_T2 + (a_1 - x_1 - (1/4)*y_1)*hbar_O2_T2;
delta_H_c_1 = abs(HR_1 - HP_1);
delta_hbar_c_1 = delta_H_c_1;
delta_h_c_1 = delta_hbar_c_1 / MW_Gas;

hc_LHV_1 = delta_h_c_1/(1/MW_Gas * (x_1*MW_CO2 + (1/2)*y_1*MW_H2O + a_1*3.76*MW_N2 + (a_1 - x_1 - (1/4)*y_1)*MW_O2));
fprintf('hc_LHV_1 = %0.2f J/kg\n', hc_LHV_1);


%============================= Question 3 =================================
% The pressure at point 3 is maximum pressure since is it pressure after
% fast combustion induced by fuel injection.

max_p = input("[Part I] Enter maximum pressure engine can sustain in MPa:")

p3 = max_p*10^6;
alpha = p3/p2;

% q_in = H_P(T3-T2) + H_P(T4-T3) - R_mix(T3-T2)
% substituting q_in = hc_LHV_1, T3 = T1*alpha*e^(n1-1), T4 = T1*alpha*beta*e^(n1-1)
% hc_LHV_1 = h_P(T3-T2) + h_P(T1*alpha*beta*e^(n1-1)-T3) - R_mix(T3-T2)
% Rearrange equation to solve for beta:

MW_products = x_1*MW_CO2 + (1/2)*y_1*MW_H2O + a_1*3.76*MW_N2 + (a_1 - x_1 - (1/4)*y_1)*MW_O2;

T3 = T1*alpha*e^(n1-1);
M_mix_1 = (x_1*MW_CO2 + (1/2)*y_1*MW_H2O + a_1*3.76*MW_N2 + (a_1 - x_1 - (1/4)*y_1)*MW_O2)/(x_1 + (1/2)*y_1 + a_1*3.76 + (a_1 - x_1 - (1/4)*y_1))
R_mix = R_u/M_mix_1;
hp_1 = abs((HP_1/MW_products)/MW_products);

beta_1 = (hc_LHV_1 - hp_1*(T3-T2) + hp_1*T3 + R_mix*(T3-T2))/(T1*hp_1*alpha*e^(n1-1));

fprintf('alpha = %0.2f\n', alpha);
fprintf('beta_1 = %0.2f\n', beta_1);


%============================= Question 4 =================================
% The formula for the net work requires the delta value.
% To find delta we need v5 and v4.
% Since this is a dual cycle, v5 = v1 and v4 = v3*beta and v3 = v2.

delta_1 = v1 / (v2*beta_1);

%for the dual cycle, the net work equation is the following
W_net_1 = 0.001*p1*v1*e^(n1-1) * (alpha*(beta_1-1) + (alpha*beta_1)/(n2-1) * (1-(1/delta_1^(n2-1))) - (1/(n1-1))*(1-(1/e^(n1-1))));

fprintf('W_net_1 = %4.2f J\n', W_net_1);


%============================= Question 5 =================================
%Find the mep of the dual cycle:
% mep = W_cycle / (v_max - v_min)
% W_cycle = W_net
% V_max = v1
% V_min = v2
mep_1 = W_net_1 / (v1 - v2);

%print mep result
fprintf('mep = %0.2f kPa\n', mep_1);


%============================= Question 6 =================================
%W_cycle, net work of cycle
%hc_LHV, lower heating value
n_th_1 = W_net_1/hc_LHV_1;
fprintf('n_th_1 = %0.2f\n', n_th_1);


%============================= Question 7 =================================
% The engine burns 1 mole of gasoline per cycle
% To generate 55000 kW or 55000 J/s:
num_cycles = 55000 / W_net_1;
num_moles = num_cycles;
mdot_1 = num_moles*0.001*MW_Gas;
fprintf('mdot_1 = %0.2f kg/s\n\n', mdot_1);


%% Part II
%============================= Question 1 =================================
% Define given values:
x_2 = 0;
y_2 = 2;

% Find all values for the chemical equation: Cx Hy   +  a(O2+ 3.76N) --> bCO2   +   cH2O  + dO2  + 3.76 aN2
a_stoich_2 = x_2 + y_2/4;
a_2 = a_stoich_2 / ER;
b_2 = x_2;
c_2 = y_2/2;
d_2 = a_2-x_2-(y_2/4);
fprintf('Part II Results:\n')
fprintf('C8H18 + %0.2f*(O2+ 3.76N) --> %0.2f*CO2 + %0.2f*H2O + %0.2f*O2 + %0.2f*N2\n', a_2, b_2, c_2, d_2, 3.76*a_2);


%============================= Question 2 =================================
% Find enthalpy of each species at T2:
% T2 is already known from part I
% The only enthalpy value which is unknown is hbar_H2_T2
MW_H2 = 2;
a1_H2_T2 = 0.03298124e02;
a2_H2_T2 = 0.08249441e-02;
a3_H2_T2 = -0.08143015e-05;
a4_H2_T2 = -0.09475434e-09;
a5_H2_T2 = 0.04134872e-11;
a6_H2_T2 = -0.10125209e04;
hbar_H2_T2 = R_u*(a6_H2_T2 + a1_H2_T2*T2^2 + (1/3)*a3_H2_T2*T2^3 + (1/4)*a4_H2_T2*T2^4 + (1/5)*a5_H2_T2*T2^5);

HR_2 = hbar_H2_T2 + a_2*hbar_O2_T2 + a_2*3.76*hbar_N2_T2;
HP_2 = (1/2)*y_2*hbar_H2O_T2 + a_2*3.76*hbar_N2_T2 + (a_2 - (1/4)*y_2)*hbar_O2_T2;

delta_H_c_2 = abs(HR_2 - HP_2);
delta_hbar_c_2 = delta_H_c_2;   % One mole of fuel is used
delta_h_c_2 = delta_hbar_c_2 / MW_H2;
hc_LHV_2 = delta_h_c_2/(1/MW_H2 * ((1/2)*y_2*MW_H2O + a_2*3.76*MW_N2 + (a_2 - (1/4)*y_2)*MW_O2));
fprintf('hc_LHV_2 = %0.2f J/kg\n', hc_LHV_2);


%============================= Question 3 =================================

% q_in = H_P(T3-T2) + H_P(T4-T3) - R_mix(T3-T2)
% substituting q_in = hc_LHV_1, T3 = T1*alpha*e^(n1-1), T4 = T1*alpha*beta*e^(n1-1)
% hc_LHV_1 = H_P(T3-T2) + H_P(T1*alpha*beta*e^(n1-1)-T3) - R_mix(T3-T2)
% Rearrange equation to solve for beta:

MW_products = x_2*MW_CO2 + (1/2)*y_2*MW_H2O + a_2*3.76*MW_N2 + (a_2 - x_2 - (1/4)*y_2)*MW_O2;

T3 = T1*alpha*e^(n1-1);
M_mix_2 = (x_2*MW_CO2 + (1/2)*y_2*MW_H2O + a_2*3.76*MW_N2 + (a_2 - x_2 - (1/4)*y_2)*MW_O2)/(x_2 + (1/2)*y_2 + a_2*3.76 + (a_2 - x_2 - (1/4)*y_2));
R_mix = R_u/M_mix_2;
hp_2 = abs((HP_2/MW_products)/MW_products);
beta_2 = (hc_LHV_2 - hp_2*(T3-T2) + hp_2*T3 + R_mix*(T3-T2))/(T1*hp_2*alpha*e^(n1-1));

fprintf('alpha = %0.2f\n', alpha);
fprintf('beta_2 = %0.2f\n', beta_2);


%============================= Question 4 =================================
% The formula for the net work requires the delta value.
% To find delta we need v5 and v4.
% Since this is a dual cycle, v5 = v1 and v4 = v3*beta and v3 = v2.

delta_2 = v1 / (v2*beta_2);

%for the dual cycle, the net work equation is the following
W_net_2 = 0.001*p1*v1*e^(n1-1) * (alpha*(beta_2-1) + (alpha*beta_2)/(n2-1) * (1-(1/delta_2^(n2-1))) - (1/(n1-1))*(1-(1/e^(n1-1))));

fprintf('W_net_2 = %4.2f J\n', W_net_2);

%============================= Question 5 =================================
%Find the mep of the dual cycle:
% mep = W_cycle / (v_max - v_min)
% W_cycle = W_net
% V_max = v1
% V_min = v2
mep_2 = W_net_2 / (v1 - v2);

%print mep result
fprintf('mep_2 = %0.2f kPa\n', mep_2);


%============================= Question 6 =================================
%W_cycle, net work of cycle
%hc_LHV, lower heating value
n_th_2 = W_net_2/hc_LHV_2;
fprintf('n_th_2 = %0.2f\n', n_th_2);


%============================= Question 7 =================================
% The engine burns 1 mole of gasoline per cycle
% To generate 55000 kW or 55000 J/s:
num_cycles = 55000 / W_net_2;
num_moles = num_cycles;
mdot_2 = num_moles*0.001*MW_H2;
fprintf('mdot_2 = %0.2f kg/s\n\n', mdot_2);
