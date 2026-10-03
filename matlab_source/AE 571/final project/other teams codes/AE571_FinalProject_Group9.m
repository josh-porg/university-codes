% AE 571 Final Project Group 9

clear;
clc;

% Part I & II

% User input required for fuel selection
fuel = input('Fuel (Gasoline/Hydrogen): ', 's');
manualInput = input('Manual Input? (Yes/No): ','s');
disp(' '); % Spacer for more consistent output formatting

% Debugging options to skip user input
%fuel = 'Gasoline';
%manualInput = 'No';

% Fuel Properties
if strcmpi(fuel, 'Gasoline')
    x = 8.26;
    y = 15.5;
else
    x = 0;
    y = 2;
end

% Parameters
if strcmpi(manualInput, 'Yes')
    % Users can enter their own parameters
    epsilon = input('Compression Ratio, ε: ');
    equivalenceRatio = input('Equivalence Ratio Ratio, ϕ: ');
    T1 = input('Ambient Temperature, T1 (K): ');
    p1 = input('Ambient Pressure, P1 (kPa): ');
    n1 = input('Polytropic Index for Compression, n1 ');
    n2 = input('Polytropic Index for Expansion, n2 ');
    p3 = input('Maximum Pressure, P3 (kPa): ');
    P = input('Power Output, P (W): ');
else
    % Paremeters included in the assignment
    epsilon = 16;
    equivalenceRatio = .8;
    T1 = 298;
    p1 = 100;
    n1 = 1.38;
    n2 = 1.25;
    p3 = 5.5 * 10^3;
    P = 55;
end


% Chemical Reaction Equation
a = (x + y/4) / equivalenceRatio;
disp('Chemical Reaction Equation: ');
if strcmpi(fuel, 'Gasoline')
    disp(cell2mat(strcat('C', num2str(x), 'H', num2str(y), {' + '}, num2str(a), {'O2 + '}, num2str(3.76*a), {'N2 → '}, num2str(x), {'CO2 + '}, num2str(y/2), {'H2O + '}, num2str(a - x - y/4), {'O2 + '}, num2str(3.76*a), 'N2')));
else
    disp(cell2mat(strcat('H', num2str(y), {' + '}, num2str(a), {'O2 + '}, num2str(3.76*a), {'N2 → '}, num2str(x), {'CO2 + '}, num2str(y/2), {'H2O + '}, num2str(a - x - y/4), {'O2 + '}, num2str(3.76*a), 'N2')));
end
disp(' '); % Spacer for more consistent output formatting

% Lower Heat Value
T2 = T1 * epsilon^(n1-1);

% Enthalpy Values
Tr = T2;
hFuelref = hBar(fuel, Tr);
hO2ref = hBar('O2', Tr);
hN2ref = hBar('N2', Tr);
hCO2ref = hBar('CO2', Tr);
hH2Oref = hBar('H2O', Tr);

hR = [hFuelref, hO2ref, hN2ref];
hP = [hCO2ref, hH2Oref, hO2ref, hN2ref];


%Combined Enthalpies
rMoles = [1, a, 3.76 * a];
pMoles = [x, y/2, a - x - y/4, 3.76 * a];
pMW = [12 + 2 * 16, 2 + 16, 2 * 16, 2 * 14];
HR = sum(hR .* rMoles);
HP = sum(hP .* pMoles);
Hc = HR - HP;
LHV = Hc / sum(pMoles .* pMW);
disp('Lower Heating Value, ∆hc (kJ/kg): ');
disp(LHV);


% Ratio of Pressure Increase
p2 = p1 * epsilon^n1;
alpha = p3/p2;
disp('Ratio of Pressure Increase, α: ');
disp(alpha);

% Ratio of Volume Increase
T3 = T2 * alpha;
Tr = T3;  

% cp calculation
cpCO2 = (hCO2ref - hBar('CO2', Tr)) / (T2 - T3);
cpH2O = (hH2Oref - hBar('H2O', Tr)) / (T2 - T3);
cpO2 = (hO2ref - hBar('O2', Tr)) / (T2 - T3);
cpN2 = (hN2ref - hBar('n2', Tr)) / (T2 - T3);

cpP = [cpCO2, cpH2O, cpO2, cpN2];
cpMix = sum(cpP .* pMoles) / sum(pMoles .* pMW);
Ru = 8.314;
cvMix = sum((cpP - Ru) .* pMoles) / sum(pMoles .* pMW);

q1 = cvMix * (T3-T2);
q2 = LHV - q1;
T4= q2 / cpMix + T3;

beta = T4/T3;
disp('Ratio of Volume Increase, β: ');
disp(beta);


% Net Cycle Work
MWair = 29;
R1 = Ru / MWair;
v1 = R1 * T1 / p1;
delta = epsilon / beta;
wCycle = p1*v1*epsilon^(n1-1) * (alpha*(beta-1) + ((alpha*beta)/(n2-1))*(1-(1/(delta^(n2-1)))) - (1/(n1-1))*(1-(1/(epsilon^(n1-1)))));
disp('Net Cycle Work, Wcycle (kJ/kg): ');
disp(wCycle);


% Cycle mep
v2 = v1 / epsilon;
mep = wCycle / (v1 - v2);
disp('Cycle mep, (kPa): ');
disp(mep);


% Thermal Efficiency
eta = wCycle / LHV;
disp('Cycle Thermal Efficiency, η: ');
disp(eta);


% Fuel Required
kg = P/wCycle;
disp('Fuel Required, (kg/s): ');
disp(kg);


% Functions

% a values

function [a1, a2, a3, a4, a5, a6] = moleculeData(molecule, T)
if strcmpi(molecule, 'Gasoline')
    a1 = -24.978;
    a2 = 256.63;
    a3 = -201.68;
    a4 = 64.75;
    a5 = .5808;
    a6 = -27.562;
elseif strcmpi(molecule, 'CO2')
    if T >1000
        a1 = .04453623 * 10^2;
        a2 = .03140168 * 10^-1;
        a3 = -.12784105 * 10^-5;
        a4 = .02393996 * 10^-8;
        a5 = -.16690333 * 10^-13;
        a6 = -.04896696 * 10^6;
    else
        a1 = .02275724 * 10^2;
        a2 = .09922072 * 10^-1;
        a3 = -.1040911 * 10^-4;
        a4 = .06866686 * 10^-7;
        a5 = -.02117280 * 10^-10;
        a6 = -.04837314 * 10^6;
    end
elseif strcmpi(molecule, 'Hydrogen')
    if T >1000
        a1 = .02991423 * 10^2;
        a2 = .07000644 * 10^-2;
        a3 = -.05633828 * 10^-6;
        a4 = -.09231578 * 10^-10;
        a5 = .15827519 * 10^-14;
        a6 = -.08350304 * 10^4;
    else
        a1 = .03298124 * 10^2;
        a2 = .08249441 * 10^-2;
        a3 = -.08143015 * 10^-5;
        a4 = -.09475434 * 10^-9;
        a5 = .04134872 * 10^-11;
        a6 = -.10125209 * 10^4;
    end
elseif strcmpi(molecule, 'H2O')
    if T >1000
        a1 = .02672145 * 10^2;
        a2 = .03056293 * 10^-1;
        a3 = -.08730260 * 10^-5;
        a4 = .12009964 * 10^-9;
        a5 = -.06391618 * 10^-13;
        a6 = -.02989921 * 10^6;
    else
        a1 = .03386842 * 10^2;
        a2 = .03474982 * 10^-1;
        a3 = -.06354696 * 10^-4;
        a4 = .06968581 * 10^-7;
        a5 = -.02506588 * 10^-10;
        a6 = -.03020811 * 10^6;
    end
elseif strcmpi(molecule, 'N2')
    if T >1000
        a1 = .02926640 * 10^2;
        a2 = .14879768 * 10^-2;
        a3 = -.05684760 * 10^-5;
        a4 = .10097038 * 10^-9;
        a5 = -.06753351 * 10^-13;
        a6 = -.09227977 * 10^4;
    else
        a1 = .03298677 * 10^2;
        a2 = .14082404 * 10^-2;
        a3 = -.03963222 * 10^-4;
        a4 = .05641515 * 10^-7;
        a5 = -.02444854 * 10^-10;
        a6 = -.10208999 * 10^4;
    end
elseif strcmpi(molecule, 'O2')
    if T >1000
        a1 = .03697578 * 10^2;
        a2 = .06135197 * 10^-2;
        a3 = -.12588420 * 10^-6;
        a4 = .01775281 * 10^-9;
        a5 = -.11364354 * 10^-14;
        a6 = -.12339301 * 10^4;
    else
        a1 = .03212936 * 10^2;
        a2 = .11274864 * 10^-2;
        a3 = -.05756150 * 10^-5;
        a4 = .13138773 * 10^-8;
        a5 = -.08768554 * 10^-11;
        a6 = -.10052490 * 10^4;
    end
end
end

% h bar values
function [h] = hBar(molecule, Tr)
Ru = 8.314;
[a1, a2, a3, a4, a5, a6] = moleculeData(molecule, Tr);
if strcmpi(molecule, 'Gasoline')
    h = 4184 * (a1 * (Tr/1000) + a2/2 * (Tr/1000)^2 + a3/3 * (Tr/1000)^3 + a4/4 * (Tr/1000)^4 - a5 * (Tr/1000)^-1 + a6);
else
    h = Ru * (a1 * Tr + a2/2 * Tr^2 + a3/3 * Tr^3 + a4/4 * Tr^4 + a5/5 * Tr^5 + a6);
end
end
