clear
clc

%% Part 1 Chemical Equation

% Atom Balance for a, b, c, d

Fuel = 'Gasoline (Heavy)';  % Defines the fuel type to extract the proper coefficients from Atom Balance
ER = 0.8;                   % Efficiency Ratio, fixed as per the engine specifications

[x, y, a, b, c, d, MW_fuel] = atombalance(Fuel, ER);

% Using the outputs of the Atom Balance Function, the coefficients of the
% chemical equation can be determined and the equation can be written.

%% Part 2 Lower Heating Value at T2

p1 = 100;       % Ambient Pressure
T1 = 298;       % Ambient Temperature
eps = 16;       % Compression Ratio
n1 = 1.38;      % Polytropic Index of Compression
n2 = 1.25;      % Polytropic Index of Expansion

p2 = p1*power(eps, n1);         % Pressure After Compression
T2 = T1*power(eps, (n1 - 1));   % Temperature After Compression

h_Fuel_ref = enthalpycalc('C8.26H15.5', 298);    % These h and cp values are calculated for the
h_O2_ref = enthalpycalc('O2', 298);              % reference temperature of 298K
h_N2_ref = enthalpycalc('N2', 298);
h_CO2_ref = enthalpycalc('CO2', 298);
h_H2O_ref = enthalpycalc('H2O', 298);

cp_fuel = (cpcalc('C8.26H15.5', T2) + cpcalc('C8.26H15.5', 298))/2;     % This section calculates the average cp
cp_O2 = (cpcalc('O2', T2) + cpcalc('O2', 298))/2;                       % for each component in the mixture, accounting
cp_N2 = (cpcalc('N2', T2) + cpcalc('N2', 298))/2;                       % for variable cp through the temperature range
cp_CO2 = (cpcalc('CO2', T2) + cpcalc('CO2', 298))/2;
cp_H2O = (cpcalc('H2O', T2) + cpcalc('H2O', 298))/2;

Hr = (h_Fuel_ref + a*h_O2_ref + 3.76*a*h_N2_ref) + (cp_fuel + a*cp_O2 + 3.76*a*cp_N2)*(T2-298);
Hp = (b*h_CO2_ref + c*h_H2O_ref + d*h_O2_ref + 3.76*a*h_N2_ref) + (b*cp_CO2 + c*cp_H2O + d*cp_O2 + 3.76*a*cp_N2)*(T2-298);

% The above lines calculate the enthalpy of the products and reactants for
% use in calculating the enthalpy change.

delta_Hc = Hr - Hp;         % Calculation of the enthalpy change between the two states

delta_hbar_c = delta_Hc;    % Since only 1 mole of fuel is present, these two quantities are equal

delta_h_c = (delta_hbar_c/(x*12 + y*1)) / ((b*44 + c*18 + d*32 + 3.76*a*28)/(x*12 + y*1)); % Calculates the LHV of the mixture per kg of mixture

%% Part 3 Ratio of Pressure Increase Alpha and Ratio of Volume Increase Beta

p3 = input('Input the Maximum Cycle Pressure in kPa, Ex: 5500: ');      % Maximum cycle pressure in kilopascals
alpha = p3/p2;                                                          % Ratio of Pressure Increase

cp_mix = (cp_fuel + a*cp_O2 + 3.76*a*cp_N2)/(MW_fuel + a*32 + 3.76*a*28);

beta = delta_h_c/(cp_mix*T2);
delta = eps/beta;

%% Part 4 Calculation of Cycle Work

Wcycle = p1*v1*(power(eps, n1 - 1))*(alpha*(beta - 1) + ((alpha*beta)/(n2 - 1))*(1 - (1/(power(delta, n2 - 1))) - (1/(n1 - 1))*(1 - (1/(power(eps, n1 - 1))))));

% The above equation calculates the total cycle work using the
% thermodynamic properties of the combustion cycle.

%% Part 5 Mean Effective Pressure of Cycle

mep = Wcycle/(v1 - v2);

% The above equation calculates the mean effective pressure using the
% simplified equation relating the cycle work and change in volume.

%% Part 6 Thermal Efficiency of the Cycle

heta = Wcycle/delta_h_c;

% The above equation calculates the thermal efficiency of the cycle using
% the simplified equation relating the cycle work to the heat input. Heat
% input is equivalent to the heat of combustion, which is why that is used
% in place of the heat input in this equation.

%% Part 7 Required Fuel Flow Rate for 55 Kilowatts Power Output

pwr = 55;
% Work in Progress

%% Atom Balance

% This function performs the balancing of the chemical equation and solves
% for the coefficients required for the equation.

function [x, y, a, b, c, d, MW_fuel] = atombalance(Fuel, ER)

if strcmpi(Fuel, 'Gasoline (Heavy)') == 1
    x = 8.26;
    y = 15.5;

    astoich = x + y/4;
    a = astoich/ER;
    b = x;
    c = y/2;
    d = a - b - (c/2);

    MW_fuel = 114.8;

elseif strcmpi(Fuel, 'Gasoline (Light)') == 1
    x = 7.76;
    y = 13.1;

    astoich = x + y/4;
    a = astoich/ER;
    b = x;
    c = y/2;
    d = a - b - (c/2);

    MW_fuel = 106.4;

elseif strcmpi(Fuel, 'Hydrogen') == 1
   x = 0;
   y = 2;

   astoich = x + y/4;
   a = astoich/ER;
   b = x;
   c = y/2;
   d = a - b - (c/2);

   MW_fuel = 2;

end

end

%% Enthalpy Calc

% This function calculates the enthalpy of each species at the required
% temperature using the specified equations.

function h = enthalpycalc(Species, T)

if strcmpi(Species, 'C8.26H15.5') == 1
    a1 = -24.078;
    a2 = 256.63;
    a3 = -201.68;
    a4 = 64.75;
    a5 = 0.5808;
    a6 = -27.562;

    theta = T/1000; 
    h = 4184*(a1*theta + a2*theta^2/2 + a3*theta^3/3 + a4*theta^4/4 + a5*theta^-2);

elseif strcmpi(Species, 'C7.76H13.1') == 1
    a1 = -22.501;
    a2 = 227.99;
    a3 = -177.26;
    a4 = 56.048;
    a5 = 0.4845;
    a6 = -17.578;

    theta = T/1000;
    h = 4184*(a1*theta + a2*theta^2/2 + a3*theta^3/3 + a4*theta^4/4 - a5*theta^-1 + a6);

elseif strcmpi(Species, 'CO2') == 1
    if T >= 1000
        a1 = 0.04453623E2;
        a2 = 0.03140168E-1;
        a3 = -0.12784105E-5;
        a4 = 0.02393996E-8;
        a5 = -0.16690333E-13;
        a6 = -0.04896696E6;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);

    elseif T < 1000
        a1 = 0.02275724E2;
        a2 = 0.09922072E-1;
        a3 = -0.10409113E-4;
        a4 = 0.06866686E-7;
        a5 = -0.0211728E-10;
        a6 = -0.04837314E6;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);
    end

elseif strcmpi(Species, 'H2O') == 1
    if T >= 1000
        a1 = 0.02672145E2;
        a2 = 0.03056293E-1;
        a3 = -0.0873026E-5;
        a4 = 0.12009964E-9;
        a5 = -0.06391618E-13;
        a6 = -0.02989921E6;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);

    elseif T < 1000
        a1 = 0.03386842E2;
        a2 = 0.03474982E-1;
        a3 = -0.06354696E-4;
        a4 = 0.06968581E-7;
        a5 = -0.02506588E-10;
        a6 = -0.03020811E6;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);
    end

elseif strcmpi(Species, 'O2') == 1
    if T >= 1000
        a1 = 0.03697578E2;
        a2 = 0.06135197E-2;
        a3 = -0.12588420E-6;
        a4 = 0.01775281E-9;
        a5 = -0.11364354E-14;
        a6 = -0.12339301E4;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);

    elseif T < 1000
        a1 = 0.03212936E2;
        a2 = 0.11274864E-2;
        a3 = -0.0575615E-5;
        a4 = 0.13138773E-8;
        a5 = -0.08768554E-11;
        a6 = -0.1005249E4;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);
    end

elseif strcmpi(Species, 'N2') == 1
    if T >= 1000
        a1 = 0.0292664E2;
        a2 = 0.14879768E-2;
        a3 = -0.0568476E-5;
        a4 = 0.10097038E-9;
        a5 = -0.06753351E-13;
        a6 = -0.09227977E4;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);

    elseif T < 1000
        a1 = 0.03298677E2;
        a2 = 0.14082404E-2;
        a3 = -0.03963222E-4;
        a4 = 0.05641515E-7;
        a5 = -0.02444854E-10;
        a6 = -0.10208999E4;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);
    end

elseif strcmpi(Species, 'H2') == 1
    if T >= 1000
        a1 = 0.02991423E2;
        a2 = 0.07000644E-2;
        a3 = -0.05633828E-6;
        a4 = -0.09231578E-10;
        a5 = 0.15827519E-14;
        a6 = -0.0835034E4;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);

    elseif T < 1000
        a1 = 0.03298124E2;
        a2 = 0.08249441E-2;
        a3 = -0.08143015E-5;
        a4 = -0.09475434E-9;
        a5 = 0.04134872E-11;
        a6 = -0.10125209E4;

        h = 8.314*(a1*T + a2/2*T^2 + a3/3*T^3 + a4/4*T^4 + a5/5*T + a6);
    end

end

end

%% CP Calc

% This function calculates the cp of each species at the specified
% temperature using the equations required for that species.

function cp = cpcalc(Species, T)

if strcmpi(Species, 'C8.26H15.5') == 1
    a1 = -24.078;
    a2 = 256.63;
    a3 = -201.68;
    a4 = 64.75;
    a5 = 0.5808;
    a6 = -27.562;

    theta = T/1000; 
    cp = 4.184*(a1 + a2*theta + a3*theta^2 + a4*theta^3 - a5*theta^-1);

elseif strcmpi(Species, 'C7.76H13.1') == 1
    a1 = -22.501;
    a2 = 227.99;
    a3 = -177.26;
    a4 = 56.048;
    a5 = 0.4845;
    a6 = -17.578;

    theta = T/1000;
    cp = 4.184*(a1 + a2*theta + a3*theta^2 + a4*theta^3 - a5*theta^-1);

elseif strcmpi(Species, 'CO2') == 1
    if T >= 1000
        a1 = 0.04453623E2;
        a2 = 0.03140168E-1;
        a3 = -0.12784105E-5;
        a4 = 0.02393996E-8;
        a5 = -0.16690333E-13;
        a6 = -0.04896696E6;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);

    elseif T < 1000
        a1 = 0.02275724E2;
        a2 = 0.09922072E-1;
        a3 = -0.10409113E-4;
        a4 = 0.06866686E-7;
        a5 = -0.0211728E-10;
        a6 = -0.04837314E6;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);
    end

elseif strcmpi(Species, 'H2O') == 1
    if T >= 1000
        a1 = 0.02672145E2;
        a2 = 0.03056293E-1;
        a3 = -0.0873026E-5;
        a4 = 0.12009964E-9;
        a5 = -0.06391618E-13;
        a6 = -0.02989921E6;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);

    elseif T < 1000
        a1 = 0.03386842E2;
        a2 = 0.03474982E-1;
        a3 = -0.06354696E-4;
        a4 = 0.06968581E-7;
        a5 = -0.02506588E-10;
        a6 = -0.03020811E6;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);
    end

elseif strcmpi(Species, 'O2') == 1
    if T >= 1000
        a1 = 0.03697578E2;
        a2 = 0.06135197E-2;
        a3 = -0.12588420E-6;
        a4 = 0.01775281E-9;
        a5 = -0.11364354E-14;
        a6 = -0.12339301E4;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);

    elseif T < 1000
        a1 = 0.03212936E2;
        a2 = 0.11274864E-2;
        a3 = -0.0575615E-5;
        a4 = 0.13138773E-8;
        a5 = -0.08768554E-11;
        a6 = -0.1005249E4;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);
    end

elseif strcmpi(Species, 'N2') == 1
    if T >= 1000
        a1 = 0.0292664E2;
        a2 = 0.14879768E-2;
        a3 = -0.0568476E-5;
        a4 = 0.10097038E-9;
        a5 = -0.06753351E-13;
        a6 = -0.09227977E4;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);

    elseif T < 1000
        a1 = 0.03298677E2;
        a2 = 0.14082404E-2;
        a3 = -0.03963222E-4;
        a4 = 0.05641515E-7;
        a5 = -0.02444854E-10;
        a6 = -0.10208999E4;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);
    end

elseif strcmpi(Species, 'H2') == 1
    if T >= 1000
        a1 = 0.02991423E2;
        a2 = 0.07000644E-2;
        a3 = -0.05633828E-6;
        a4 = -0.09231578E-10;
        a5 = 0.15827519E-14;
        a6 = -0.0835034E4;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);

    elseif T < 1000
        a1 = 0.03298124E2;
        a2 = 0.08249441E-2;
        a3 = -0.08143015E-5;
        a4 = -0.09475434E-9;
        a5 = 0.04134872E-11;
        a6 = -0.10125209E4;

        cp = 8.314*(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4);
    end

end

end