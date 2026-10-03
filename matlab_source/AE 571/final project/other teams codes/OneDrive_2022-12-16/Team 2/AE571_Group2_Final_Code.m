%% PART 1, Gasoline as a Fuel

clc
clear

fprintf('<strong>Gasoline: </strong>\n')

    % Givens
        epsilon = 16;
        ER = 0.8;
        T1 = 298;
        p1 = 100;
        n1 = 1.38;
        n2 = 1.25;
        Ru = 8.314; % kJ/kmol K

% 1) Determine Chemical Reaction Equation

        [a,x,y] = atombalance('Gasoline',ER); % Calling the function for the atom balance
        fprintf("C8.26H15.5 + %3.2f*(O2 + 3.76N2) --> %3.2fCO2 + %3.2fH2O(g) + %3.2fO2 + %3.2f*3.76N2 \n",a,x,y/2,(a-x-y/4),a)

    % Defining Mixture Molecular Weights
        [MW_fuel,MW_mix_r,MW_mix_p] = mixMW(a,x,y);

    % Calculate Properties at Point 2
        p2 = p1*epsilon^(n1);   % Polytropic Relationship
        T2 = T1*epsilon^(n1-1); % Ideal Gas Law

% 2) The lower heat value of the fuel at T2 (the temperature at the end of compression) 
    % Call a function to calculate the total enthalpy for either reactants or products
        Hr = mixenth('R','Gasoline',T2,a,x,y); % kJ/kmol
        Hp = mixenth('P','Gasoline',T2,a,x,y); % kJ/kmol

    % Finding the LHV in the correct units
        dHc = Hr-Hp; % kJ/kmol
        dhbar = dHc/1; % kJ/kmol fuel
        dhc = dhbar/MW_fuel; % kJ/kg fuel
        LHV = dhc/((MW_mix_p*(x+(y/2)+(a-x-(y/4))+a*3.76))/MW_fuel); %kJ/kg total
        fprintf('The lower heating value of the cycle is %3.2f kJ/kg\n', LHV)

% 3) Alpha and beta if the maximum pressure the engine can sustain is 5.5 MPa 

    % Calculating properties p3, p4, and T3
        p3 = 5.5*10^3; % Max pressure occurs at point 3
        p4 = p3; % This is true because we are doing a dual cycle
        alpha = p3/p2;
        T3 = T2*alpha; % Using ideal gas law

    % solving for beta
        nmols_p = x+(y/2)+(a-x-(y/4))+3.76*a; % Number of moles of products 
        Rmix = Ru/MW_mix_p; % kJ/kg*K
        dummy = LHV + mixenth('P',0,T2,a,x,y)/(MW_mix_p*nmols_p) + Rmix*(T3-T2); % Dummy variable, equal to enthalpy of products at T4 used for the loop
        T4 = T3; % Set this to start the loop, then iterate to find what T4 actually is
        h4 = mixenth('P',0,T4,a,x,y)/(MW_mix_p*nmols_p); %kJ

        while abs(dummy-h4)>= 5
            if dummy > h4
                T4 = T4 + 0.01; % Increase T4 if needed
            else
                T4 = T4 - 0.01; % Or decrease T4 if needed
            end
            h4 = mixenth('P',0,T4,a,x,y)/(MW_mix_p*nmols_p); % Enthalpy at point 4
        end

        beta = T4/T3; % T4/T3 = V4/V3 from the ideal gas law since P3 = P4
        R_reactants = Ru/MW_mix_r;
        V1 = R_reactants*T1/p1; % Ideal gas law
        delta = epsilon/beta; % delta = V5/V4
        fprintf('The pressure increase of the cycle, alpha, is %3.3f\n', alpha)
        fprintf('The volume increase, beta, of the cycle is %3.3f\n', beta)

% 4) Net work of the cycle
        Wcycle = p1*V1*(epsilon^(n1-1))*(alpha*(beta-1)+(((alpha*beta)/(n2-1)))*(1-(1/delta^(n2-1)))-(1/(n1-1))*(1-(1/epsilon^(n1-1))));
        fprintf('The net work of the cycle is %3.2f kJ\n', Wcycle)

% 5) Mean effective pressure of the cycle
        mep = Wcycle/(V1*(1-(1/epsilon)));
        fprintf('The mean effective pressure of the cycle is %3.2f kPa\n', mep)

% 6) Thermal efficiency of the cycle
        nth = 100*Wcycle/LHV;
        fprintf('The Thermal efficieny of the cycle is %3.2f%%\n', nth)

% 7) How many kg of fuel per sec needs to be injected to prudece 55kW of power
    % kW = kJ/sec = (kJ/kgfuel)*(kgfuel/sec)
        power = 55; % kW % Given 
        fuelrate = power/dhc; % kg fuel/sec, make sure units work out 
        fprintf('It will take %3.6f kg of fuel per second to produce 55 kW of power\n\n', fuelrate)
%% Part 2, Hydrogen as a fuel (H2)

fprintf('<strong>Hydrogen: </strong>\n')

    % Givens
        epsilon = 16;
        ER = 0.8;
        T1 = 298;
        p1 = 100;
        n1 = 1.38;
        n2 = 1.25;
        Ru = 8.314; % kJ/kmol K

% 1) Determine Chemical Reaction Equation

        [a,x,y] = atombalance('Hydrogen',ER);
        fprintf("H2 + %3.2f*(O2 + 3.76N2) --> %3.2fCO2 + %3.2fH2O(g) + %3.2fO2 + %3.2f*3.76N2 \n",a,x,y/2,(a-x-y/4),a)
    
    % Defining Mixture Molecular Weights
        [MW_fuel,MW_mix_r,MW_mix_p] = mixMW(a,x,y);

    % Calculate Properties at Point 2
        p2 = p1*epsilon^(n1);   % Polytropic Relationship
        T2 = T1*epsilon^(n1-1); % Ideal gas law

% 2) The lower heat value of the fuel at T2 (the temperature at the end of compression) 
    % Call a function to calculate the total enthalpy of either reactants or products 
        Hr = mixenth('R','Hydrogen',T2,a,x,y); % kJ/kmol
        Hp = mixenth('P','Hydrogen',T2,a,x,y); % kJ/kmol

    % Finding the LHV in the correct units 
        dHc = Hr-Hp; % kJ/kmol
        dhbar = dHc/1; % kJ/kmol fuel
        dhc = dhbar/MW_fuel; % kJ/kg fuel
        LHV = dhc/((MW_mix_p*(x+(y/2)+(a-x-(y/4))+a*3.76))/MW_fuel); % kJ/kg total
        fprintf('The lower heating value of the cycle is %3.2f kJ/kg\n', LHV)

% 3) Alpha and beta if the maximum pressure the engine can sustain is 5.5 MPa 

    % Calculating properties p3, p4, and T3
        p3 = 5.5*10^3; % Max pressure occurs at point 3 
        p4 = p3; % This is true because we have a dual cycle
        alpha = p3/p2;
        T3 = T2*alpha; % Ideal gas law

    % solve for beta
        nmols_p = x+(y/2)+(a-x-(y/4))+3.76*a; % Number of moles of products 
        Rmix = Ru/MW_mix_p; % kJ/kg K
        dummy = LHV + mixenth('P',0,T2,a,x,y)/(MW_mix_p*nmols_p) + Rmix*(T3-T2); % Dummy variable, equal to enthalpy of products at T4 for the loop
        T4 = T3; % Set this equal to start the loop, then iterate to find T4
        h4 = mixenth('P',0,T4,a,x,y)/(MW_mix_p*nmols_p); % kJ

        while abs(dummy-h4)>= 5
            if dummy > h4
                T4 = T4 + 0.01; % Increase T4 if needed
            else
                T4 = T4 - 0.01; % Decrease T4 if needed
            end
            h4 = mixenth('P',0,T4,a,x,y)/(MW_mix_p*nmols_p); % Enthlapy at point 4
        end

        beta = T4/T3; % T4/T3 = V4/V3 also from the ideal gas law
        R_reactants = Ru/MW_mix_r;
        V1 = R_reactants*T1/p1; % Ideal gas law
        delta = epsilon/beta; % delta = V5/V4
        fprintf('The pressure increase of the cycle, alpha, is %3.3f\n', alpha)
        fprintf('The volume increase, beta, of the cycle is %3.3f\n', beta)

% 4) Net work of the cycle
        Wcycle = p1*V1*(epsilon^(n1-1))*(alpha*(beta-1)+(((alpha*beta)/(n2-1)))*(1-(1/delta^(n2-1)))-(1/(n1-1))*(1-(1/epsilon^(n1-1))));
        fprintf('The net work of the cycle is %3.2f kJ\n', Wcycle)

% 5) Mean effective pressure of the cycle 
        mep = Wcycle/(V1*(1-(1/epsilon)));
        fprintf('The mean effective pressure of the cycle is %3.2f kPa\n', mep)

% 6) Thermal efficiency of the cycle
        nth = 100*Wcycle/LHV; % Multiplies by 100 the get the percentage instead of the decimal
        fprintf('The Thermal efficieny of the cycle is %3.2f%%\n', nth)

% 7) How many kg of fuel per sec needs to be injected to prudece 55kW of power
    % kW = kJ/sec = (kJ/kgfuel)*(kgfuel/sec)
        power = 55; % kW
        fuelrate = power/dhc; % kg fuel/sec, make sure units work out
        fprintf('It will take %3.6f kg of fuel per second to produce 55 kW of power\n', fuelrate)
%% Functions

function [a, x, y] = atombalance(Fuel,ER)

if strcmpi(Fuel,'Hydrogen') == 1
    x = 0;
    y = 2;
elseif strcmpi(Fuel,'Gasoline') == 1
    x = 8.26;
    y = 15.5;

end
a_stoich= x + y/4;
a = a_stoich/ER;
end

function [MW_fuel, MW_mix_r, MW_mix_p] = mixMW(a,x,y)

% Defining Molecular Weights
MW_H = 1;
MW_O = 16;
MW_C = 12;
MW_N = 14;
MW_fuel = x*MW_C + y*MW_H;

MW_mix_r = (MW_fuel+a*(2*MW_O+2*3.76*MW_N))/(a*(1+3.76)+1);
MW_mix_p = (x*(MW_C+2*MW_O)+(y/2)*(MW_H*2+MW_O)+(a-x-(y/4))*(MW_O*2)+a*MW_N*2*3.76)/(x+(y/2)+(a-x-(y/4))+a*3.76);

end

function h = enthalpies(Species, T)

Ru = 8.314; %KJ/Kmol * K
    % Fuel Species Enthalpies
if strcmpi(Species,"Gasoline") == 1
    a1 = -24.078;
    a2 = 256.63;
    a3 = -201.68;
    a4 = 64.750;
    a5 = 0.5808;
    a6 = -27.562;
    % Species Temp-dependent Enthalpies
elseif strcmpi(Species,"CO") == 1
    if T >= 1000 
        a1 = 0.03025079e+2;
        a2 = 0.14426885e-2;
        a3 = -.05630827e-5;
        a4 = 0.10185813e-9;
        a5 = -.06910951e-13;
        a6 = -.14268350e+5;
    elseif T < 1000
        a1 = 0.03262451e+2;
        a2 = 0.15119409e-2;
        a3 = -.03881755e-4;
        a4 = 0.05581944e-7;
        a5 = -.02474951e-10;
        a6 = -.14310539e+5;
    end

elseif strcmpi(Species,"CO2") == 1
    if T >= 1000
        a1 = 0.04453623e+2;
        a2 = 0.03140168e-1;
        a3 = -.12784105e-5;
        a4 = 0.02393996e-8;
        a5 = -.16690333e-13;
        a6 = -.04896696e+6;
    elseif T < 1000
        a1 = 0.02275724e+2;
        a2 = 0.09922072e-1;
        a3 = -.10409113e-4;
        a4 = 0.06866686e-7;
        a5 = -.02117180e-10;
        a6 = -.04837314e+6;
    end

elseif strcmpi(Species,"N2") == 1
    if T >= 1000
        a1 = 0.02926640e+2;
        a2 = 0.14879768e-2;
        a3 = -.05684760e-5;
        a4 = 0.10097038e-9;
        a5 = -.06753351e-13;
        a6 = -.09227977e+4;
    elseif T < 1000
        a1 = 0.03298677e+2;
        a2 = 0.14082404e-2;
        a3 = -.03963222e-4;
        a4 = 0.05641515e-7;
        a5 = -.02444854e-10;
        a6 = -.10208999e+4;
    end

elseif strcmpi(Species,"O2") == 1
    if T >= 1000
        a1 = 0.03697578e+2;
        a2 = 0.06135197e-2;
        a3 = -.12588420e-6;
        a4 = 0.01775281e-9;
        a5 = -.11364354e-14;
        a6 = -.12339301e+4;
    elseif T < 1000
        a1 = 0.03212936e+2;
        a2 = 0.11274864e-2;
        a3 = -.05756150e-5;
        a4 = 0.13138773e-8;
        a5 = -.08768554e-11;
        a6 = -.10052490e+4;
    end
elseif strcmpi(Species,"H2O") == 1
    if T >= 1000
        a1 = 0.02672145e+2;
        a2 = 0.03056293e-1;
        a3 = -.08730260e-5;
        a4 = 0.12009964e-9;
        a5 = -.06391618e-13;
        a6 = -.02989921e+6;
    elseif T < 1000
        a1 = 0.03386842e+2;
        a2 = 0.03474982e-1;
        a3 = -.06354686e-4;
        a4 = 0.06968581e-7;
        a5 = -.02506588e-10;
        a6 = -.03020811e+6;
    end
elseif strcmpi(Species,"Hydrogen") == 1
    if T >= 1000
        a1 = 0.02991423e+2;
        a2 = 0.07000644e-2;
        a3 = -.05633828e-6;
        a4 = -.09231578e-10;
        a5 = 0.15827519e-14;
        a6 = -.08350340e+4;
    elseif T < 1000
        a1 = 0.03298124e+2;
        a2 = 0.08249441e-2;
        a3 = -.08143015e-5;
        a4 = -.09475434e-9;
        a5 = 0.04134872e-11;
        a6 = -.10125209e+4;
    end
end

if strcmpi(Species,"Gasoline") == 1
    theta = T/1000;
    h_i = (4184*(a1*theta + a2*(theta^2)/2 + a3*(theta^3)/3 + a4*(theta^4)/4 - a5/theta + a6)); % Enthalpy of formation
    Cp = @(x) (4.184*(a1 + a2*x/1000 + a3*(x/1000).^2 + a4*(x/1000).^3 + a5*(x/1000).^-2));
    Cpdt = integral(Cp, 298, T); % Integrate to get specific enthalpy

    h = (h_i + Cpdt); % kJ/kg*kmol
else
    h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 +(a5/5)*T^5);
end
end

function H = mixenth(RorP,Fuel,T,a,x,y)
    if strcmpi(RorP,'R') == 1
        H = enthalpies(Fuel,T) + a*(enthalpies('O2',T)+3.76*enthalpies('N2',T)); %kJ
    elseif strcmpi(RorP,'P') == 1
        H = x*enthalpies('CO2',T)+(y/2)*enthalpies('H2O',T)+(a-x-(y/4))*enthalpies('O2',T)+a*3.76*enthalpies('N2',T); %kJ
    else 
        return
    end
end