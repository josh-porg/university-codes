% AE 571 Final Project Code
% Team Name: The Super Big Canards
% Lucas Russell, Dean Foster, Rhett Heide, Jeb Marshall, Camden Caulfield,
% Charlie Platt, Kadin Olson, Niels Braaten


% ** TO RUN IF WANTING TO OBTAIN SOLUTIONS FOR THE ORIGINAL ASSIGNED
% PROBLEM:
% TASK 1: Press run. Type no into command window, press enter. Type gasoline into
% command window, press enter. Type 5.5 into command window, press enter.
% TASK 2: Press run. Type no into command window, press enter. Type H2 into
% command window, press enter. Type 5.5 into command window, press enter.



% Have a good day, my friend.

%% Task 1,2 Part 1
% Determines the reaction equation and loads in needed information about
% the species and engine being analyzed. Produces a balanced chemical
% equation of the combustion process to be used in analysis.

clear, clc, close

tref = 298; %reference temp, K
Ru = 8.31451; % universal gas constant
thetaref = tref/1000;
er = 0.8; % equivalence ratio
epsilon = 16; % expansion ratio
nfuel = 1; % moles of fuel
p1 = 100000; % pressure at state 1
n1 = 1.38; % polytropic compression index
n2 = 1.25; % polytropic expansion index
T = tref*epsilon^(n1-1); % temperature at state 2
cfuel = input('Custom fuel? (input yes or no, fuels included in this program are methane, propane, gasoline, H2, and diesel): ', 's'); % user input fuel
MW_o2 = 31.9988; % molecular weights and gas constant of air
MW_co2 = 44.01;
MW_n2 = 28.0134;
MW_h2o = 20.1588;
R_air = 287;

theta = T/1000;
% specific heat polynomial coefficients for 300-1000 kelvin
a1o2r = 0.03212936E2;
a2o2r = 0.11274864E-2;
a3o2r = -0.0575615E-5;
a4o2r = 0.13138773E-8;
a5o2r = -0.08768554E-11;
a6o2r = -0.1005249E4;

a1co2r = 0.02275724E2;
a2co2r = 0.09922072E-1;
a3co2r = -0.10409113E-4;
a4co2r = 0.06866686E-7;
a5co2r = -0.0211728E-10;
a6co2r = -0.04837314E6;

a1h2or = 0.03386842E2;
a2h2or = 0.03474982E-1;
a3h2or = -0.06354696E-4;
a4h2or = 0.06968581E-7;
a5h2or = -0.02506588E-10;
a6h2or = -0.03020811E6;

a1n2r = 0.03298677E2;
a2n2r = 0.14082404E-2;
a3n2r = -0.03963222E-4;
a4n2r = 0.05641515E-7;
a5n2r = -0.02444854E-10;
a6n2r = -0.10208999E4;

% gas coefficients at reference temp

% Fuel property selection dialog. Changes fuel properties according to user input.
if strcmpi('No',cfuel) == 1
    fuel = input('Fuel Name (H2, gasoline, methane, propane, diesel): ','s');
    if strcmpi('methane',fuel) == 1
        MW_f = 16.043;
        x = 1;
        y = 4;
        a1f = -0.29149;
        a2f = 26.327;
        a3f = -10.610;
        a4f = 1.5656;
        a5f = 0.16573;
        a6f = -18.331;
        cpfuel = 4.184*(a1f + a2f*theta + a3f*theta^2 + a4f*theta^3 + a5f*theta^-2);
        hfuel_i = 4184*(a1f*theta + a2f*(theta^2/2) + a3f*(theta^3/3) + a4f*(theta^4/4) - a5f*(theta^-1) + a6f);
    elseif strcmpi('propane',fuel) == 1
        MW_f = 44.096;
        x = 3;
        y = 8;
        a1f = -1.4867;
        a2f = 74.339;
        a3f = -39.065;
        a4f = 8.0543;
        a5f = 0.01219;
        a6f = -27.313;
        cpfuel = 4.184*(a1f + a2f*theta + a3f*theta^2 + a4f*theta^3 + a5f*theta^-2);
        hfuel_i = 4184*(a1f*theta + a2f*(theta^2/2) + a3f*(theta^3/3) + a4f*(theta^4/4) - a5f*(theta^-1) + a6f);
    elseif strcmpi('diesel',fuel) == 1
        MW_f = 148.6;
        x = 10.8;
        y = 18.7;
        a1f = -9.1063;
        a2f = 246.97;
        a3f = -143.74;
        a4f = 32.329;
        a5f = 0.0518;
        a6f = -50.128;
        cpfuel = 4.184*(a1f + a2f*theta + a3f*theta^2 + a4f*theta^3 + a5f*theta^-2);
        hfuel_i = 4184*(a1f*theta + a2f*(theta^2/2) + a3f*(theta^3/3) + a4f*(theta^4/4) - a5f*(theta^-1) + a6f);
    elseif strcmpi('gasoline', fuel) == 1
        MW_f = 114.8;
        x = 8.26;
        y = 15.5;
        a1f = -24.078;
        a2f = 256.63;
        a3f = -201.68;
        a4f = 64.75;
        a5f = 0.5808;
        a6f = -27.562;
        cpfuel = 4.184*(a1f + a2f*theta + a3f*theta^2 + a4f*theta^3 + a5f*theta^-2);
        hfuel_i = 4184*(a1f*theta + a2f*(theta^2/2) + a3f*(theta^3/3) + a4f*(theta^4/4) - a5f*(theta^-1) + a6f);
    elseif strcmpi('H2',fuel) == 1
        MW_f = 2;
        x = 0;
        y = 2;
        a1h2r = 0.03298124E2;
        a2h2r = 0.08249441E-2;
        a3h2r = -0.08143015E-5;
        a4h2r = -0.09475434E-9;
        a5h2r = 0.04134872E-11;
        a6h2r = -0.10125209E4;
        hfuel_i = Ru*(a6h2r + a1h2r*tref + a2h2r*(tref^2/2) + a3h2r*(tref^3/3) + a4h2r*(tref^4/4) + a5h2r*(tref^5/5));
    end
elseif strcmpi('Yes',cfuel) == 1
x = input('Number of Carbon Atoms in Hydrocarbon Fuel: ');
y = input('Number of Hydrogen Atoms in Hydrocarbon Fuel: ');
a1f = input('a1 coefficient: ');
a2f = input('a2 coefficient: ');
a3f = input('a3 coefficient: ');
a4f = input('a4 coefficient: ');
a5f = input('a5 coefficient: ');
a6f = input('a6 coefficient: ');
end

% polynomial coefficients for T> 1000 kelvin
if T >= 1000
a1o2 = 0.03697578E2;
a2o2 = 0.06135197E-2;
a3o2 = -0.1258842E-6;
a4o2 = 0.01775281E-9;
a5o2 = -0.11364354E-14;
a6o2 = -0.12339301E4;



a1co2 = 0.04453623E2;
a2co2 = 0.03140168E-1;
a3co2 = -0.12784105E-5;
a4co2 = 0.02393996E-8;
a5co2 = -0.16690333E-13;
a6co2 = -0.04896696E6;



a1h2o = 0.02672145E2;
a2h2o = 0.03056293E-1;
a3h2o = -0.0873026E-5;
a4h2o = 0.12009964E-9;
a5h2o = -0.06391618E-13;
a6h2o = -0.02989921E6;



a1n2 = 0.0292664E2;
a2n2 = 0.14879768E-2;
a3n2 = -0.0568476E-5;
a4n2 = 0.10097038E-9;
a5n2 = -0.06753351E-13;
a6n2 = -0.09227977E4;

a1h2 = 0.02991423E2;
a2h2 = 0.07000644E-2;
a3h2 = -0.05633828E-06;
a4h2 = -0.09231578E-10;
a5h2 = 0.15827519E-14;
a6h2 = -0.08350340E4;

% polynomial coefficients for T < 1000 kelvin
elseif T < 1000
a1o2 = 0.03212936E2;
a2o2 = 0.11274864E-2;
a3o2 = -0.0575615E-5;
a4o2 = 0.13138773E-8;
a5o2 = -0.08768554E-11;
a6o2 = -0.1005249E4;

a1co2 = 0.02275724E2;
a2co2 = 0.09922072E-1;
a3co2 = -0.10409113E-4;
a4co2 = 0.06866686E-7;
a5co2 = -0.0211728E-10;
a6co2 = -0.04837314E6;

a1h2o = 0.03386842E2;
a2h2o = 0.03474982E-1;
a3h2o = -0.06354696E-4;
a4h2o = 0.06968581E-7;
a5h2o = -0.02506588E-10;
a6h2o = -0.03020811E6;

a1n2 = 0.03298677E2;
a2n2 = 0.14082404E-2;
a3n2 = -0.03963222E-4;
a4n2 = 0.05641515E-7;
a5n2 = -0.02444854E-10;
a6n2 = -0.10208999E4;


a1h2 = 0.03298124E2;
a2h2 = 0.08249441E-2;
a3h2 = -0.08143015E-5;
a4h2 = -0.09475434E-9;
a5h2 = 0.04134872E-11;
a6h2 = -0.10125209E4;
end





% calculating enthalpies at state 1 (i) and state 2 (f)
ho2_f = Ru*(a6o2 + a1o2*T + a2o2*(T^2/2) + a3o2*(T^3/3) + a4o2*(T^4/4) + a5o2*(T^5/5));
ho2_i = Ru*(a6o2r + a1o2r*tref + a2o2r*(tref^2/2) + a3o2r*(tref^3/3) + a4o2r*(tref^4/4) + a5o2r*(tref^5/5));

hco2_i = Ru*(a6co2 + a1co2*tref + a2co2*(tref^2/2) + a3co2*(tref^3/3) + a4co2*(tref^4/4) + a5co2*(tref^5/5));
hco2_f = Ru*(a6co2 + a1co2*T + a2co2*(T^2/2) + a3co2*(T^3/3) + a4co2*(T^4/4) + a5co2*(T^5/5));

hh2o_i = Ru*(a6h2o + a1h2o*tref + a2h2o*(tref^2/2) + a3h2o*(tref^3/3) + a4h2o*(tref^4/4) + a5h2o*(tref^5/5));
hh2o_f = Ru*(a6h2o + a1h2o*T + a2h2o*(T^2/2) + a3h2o*(T^3/3) + a4h2o*(T^4/4) + a5h2o*(T^5/5));

hn2_f = Ru*(a6n2 + a1n2*T + a2n2*(T^2/2) + a3n2*(T^3/3) + a4n2*(T^4/4) + a5n2*(T^5/5));
hn2_i = Ru*(a6n2r + a1n2r*tref + a2n2r*(tref^2/2) + a3n2r*(tref^3/3) + a4n2r*(tref^4/4) + a5n2r*(tref^5/5));




% Finding atom balance
a_stoich = x + y/4;
a = a_stoich./er;
b = x;
c = y/2;
d = a - b - (c/2);
e = 3.76*a;

% prints the chemical reaction equation of the combustion process.
fprintf("\nThe Chemical Reaction Equation is C_%2.2fH_%2.2f + %2.2f(O_2 + 3.76N_2) ---> %2.2fCO_2 + %2.2fH2O + %2.2fO_2 + %2.2fN_2 \n", x,y,a,b,c,d,e)

%% Task 1,2 Part 2

% calculates the lower heating value of the fuel chosen in part 1, in kJ/kg
% of mixture.

% enthalpy difference between reactants and products
delta_H = (hfuel_i + a*(ho2_i+3.76*hn2_i)) - (b*hco2_i + c*hh2o_i + d*ho2_f + e*hn2_i);

delta_hbar = delta_H/nfuel; % molar specific enthalpy difference

delta_hc = delta_hbar/MW_f; % lower heating value, kJ/kg of fuel

LHV = delta_hc*(1/((1/MW_f)*(b*MW_co2 + c*MW_h2o + d*MW_o2 + e*MW_n2))) % lower heating value, kJ/kg of mixture
%% Task 1,2 Part 3

% calculates alpha and beta of the cycle using previous data and user
% inputed maximum allowable pressure.

% Maximum allowable pressure input prompt, MPa. Input 5.5 for assigned problem, but it
% can be changed to any arbitrary value.
p3 = input('Maximum allowable cycle pressure (MPa): ')*10^6;
v1 = (R_air*tref)/p1;
p2 = p1*(epsilon)^n1;
alpha = p3/p2; % alpha for a dual cycle.
T3 = alpha*T;

% more coefficients. These are needed again for the coming piecewise
% symbolic functions of cp and cv.
a1o2 = 0.03697578E2;
a2o2 = 0.06135197E-2;
a3o2 = -0.1258842E-6;
a4o2 = 0.01775281E-9;
a5o2 = -0.11364354E-14;
a6o2 = -0.12339301E4;

a1co2 = 0.04453623E2;
a2co2 = 0.03140168E-1;
a3co2 = -0.12784105E-5;
a4co2 = 0.02393996E-8;
a5co2 = -0.16690333E-13;
a6co2 = -0.04896696E6;

a1h2o = 0.02672145E2;
a2h2o = 0.03056293E-1;
a3h2o = -0.0873026E-5;
a4h2o = 0.12009964E-9;
a5h2o = -0.06391618E-13;
a6h2o = -0.02989921E6;

a1n2 = 0.0292664E2;
a2n2 = 0.14879768E-2;
a3n2 = -0.0568476E-5;
a4n2 = 0.10097038E-9;
a5n2 = -0.06753351E-13;
a6n2 = -0.09227977E4;

a1h2 = 0.02991423E2;
a2h2 = 0.07000644E-2;
a3h2 = -0.05633828E-06;
a4h2 = -0.09231578E-10;
a5h2 = 0.15827519E-14;
a6h2 = -0.08350340E4;

% creating symbolic functions of the cp of each individual species in the
% products
syms temp
hn2 = Ru*(a6n2 + a1n2*temp + a2n2*(temp^2/2) + a3n2*(temp^3/3) + a4n2*(temp^4/4) + a5n2*(temp^5/5));
hn2ref = Ru*(a6n2r + a1n2r*temp + a2n2r*(temp^2/2) + a3n2r*(temp^3/3) + a4n2r*(temp^4/4) + a5n2r*(temp^5/5));
hn2final = piecewise(temp<1000,hn2ref,temp>=1000, hn2);
cpn2 = diff(hn2final);

ho2 = Ru*(a6o2 + a1o2*temp + a2o2*(temp^2/2) + a3o2*(temp^3/3) + a4o2*(temp^4/4) + a5o2*(temp^5/5));
ho2ref = Ru*(a6o2r + a1o2r*temp + a2o2r*(temp^2/2) + a3o2r*(temp^3/3) + a4o2r*(temp^4/4) + a5o2r*(temp^5/5));
ho2final = piecewise(temp<1000,ho2ref,temp>=1000, ho2);
cpo2 = diff(ho2final);

hh2o = Ru*(a6h2o + a1h2o*temp + a2h2o*(temp^2/2) + a3h2o*(temp^3/3) + a4h2o*(temp^4/4) + a5h2o*(temp^5/5));
hh2oref = Ru*(a6h2or + a1h2or*temp + a2h2or*(temp^2/2) + a3h2or*(temp^3/3) + a4h2or*(temp^4/4) + a5h2or*(temp^5/5));
hh2ofinal = piecewise(temp<1000,hh2oref,temp>=1000, hh2o);
cph2o = diff(hh2ofinal);

hco2 = Ru*(a6co2 + a1co2*temp + a2co2*(temp^2/2) + a3co2*(temp^3/3) + a4co2*(temp^4/4) + a5co2*(temp^5/5));
hco2ref = Ru*(a6co2r + a1co2r*temp + a2co2r*(temp^2/2) + a3co2r*(temp^3/3) + a4co2r*(temp^4/4) + a5co2r*(temp^5/5));
hco2final = piecewise(temp<1000,hco2ref,temp>=1000, hco2);
cpco2 = diff(hco2final);

% combining the cp piecewise functions into one higher level function to
% represent the cp of the entire mixture with respect to temperature.
cp_mix = (b*cpco2 + c*cph2o + d*cpo2 + e*cpn2)/(b*MW_co2 +c*MW_h2o + d*MW_o2 + e*MW_n2);
% mole fractions of species in products
xco2 = b/(b+c+d+e);
xo2 = d/(b+c+d+e);
xh2o = c/(b+c+d+e);
xn2 = e/(b+c+d+e);
% molecular weight of the mixture
MW_mix = (MW_n2*xn2 + MW_o2*xo2 + MW_h2o*xh2o + MW_co2*xco2);
% gas constant of the mixture
R_mix = (b*Ru + c*Ru + d*Ru + e*Ru)/(b*MW_co2 +c*MW_h2o + d*MW_o2 + e*MW_n2);
% cv function of the mixture
cv_mix = cp_mix-R_mix;

% integrating over the cv function to calculate constant volume heat input
q23 = vpa(int(cv_mix,T,T3));
% solving algebraically for constant pressure heat input in a dual cycle.
q34 = LHV-q23;
% calculating T4 by symbolic integration and filtering out mathematical solutions
% that do not reflect reality.
syms x
T4 = vpa(solve(int(cp_mix,T3,x) == q34));
T4 = T4(2);

% property calculations
v2 = v1/epsilon;
v3 = v2;
p4 = p3;
v4 = (1000*R_mix*T4)/p4;

% alpha and beta calculations. First block is for special scenario of otto
% cycle. Second block is for a normal dual cycle. Third block is for
% special scenario of diesel cycle. Which scenario is used depends on the
% user's maximum allowable pressure input. Block 4 will tell the user to
% input a higher pressure value since it is not physical to have an
% allowable pressre below pressure at the end of compression.
if v4 <= v3 & p3 > p2
    beta = 1
    q23 = LHV;
    q34 = 0;
    T3 = vpa(solve(int(cv_mix,T,x) == q23));
    T3 = T3(1);
    alpha = T3/T
    p3 = alpha*p2;
elseif v4 > v3 & p3 > p2
    beta = v4/v3
    alpha = p3/p2
elseif p3 == p2
    q23 = 0;
    q34 = LHV;
    T3 = vpa(solve(int(cp_mix,T,x) == q34));
    T3 = T3(2);
    v3 = (1000*R_mix*T3)/p3;
    beta = v3/v2
    alpha = 1
elseif p3 < p2
    disp('Maximum Allowable Pressure is too low. Please input a higher value.')
end

%% Task 1,2 Part 4
% calculates net work of the cycle in kJ/kg
% calculating expansion ratio
delta = epsilon/beta;

% calculating the net cycle work in kJ/kg
W_cycle = (p1*v1*epsilon^(n1-1))*(alpha*(beta-1)+((alpha*beta)/(n2-1))*(1-1/(delta^(n2-1)))-((1/(n1-1))*(1-1/(epsilon^(n1-1)))))/1000

%% Task 1,2 Part 5
% calculates mean effective pressure in MPa

mep = ((W_cycle*1000)/(v1-v2))/1E6

%% Task 1,2 Part 6

Thermal_efficiency = W_cycle/LHV

%% Task 1,2 Part 7
% calculates required mass flow rate of fuel to sustain 55 kW power, in
% kg/s

% Required power, kW
P_required = 55;

% Required mass flow rate of fuel, kg/s.
mdot_fuel = P_required/delta_hc
