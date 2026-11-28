clc;
clearvars;

% AE 571 Final Group Project

% Given
CR = 16;    % compression ratio
ER = 0.8;   % equivalence ratio 
T1 = 298;   % ambient temperature (K)
P1 = 100;   % ambient pressure (kPa)
n1 = 1.38;  % polytropic index for compression
n2 = 1.25;  % polytropic index for expansion


% Part I: Gasoline used as fuel
fuel = 'C8.26H15.5';    % Fuel chosen

% Q1: Chemical reaction equation
[a,b,c,d,x,y,MW_fuel] = atombalancefinal(fuel,ER);      % calculates coefficents for chemical reaction equation and molecular weight of fuel

% Q2: The lower heating value of the fuel at T2 (temperature at the end of
% compression)
P2 = P1*(CR^n1);            % Pressure at point 2 (kPa)
T2 = T1*(CR^(n1-1));        % Temperature at point 2 (K)

% % Enthalpies
% Reactants
h_fuel = enthalpycalculatorfinal(fuel,T1);  % enthalpy of fuel at T1 (kJ/kmol)
h_O2_R = enthalpycalculator2('O2',T1);      % enthalpy of O2 at T1 (kJ/kmol)
h_N2_R = enthalpycalculator2('N2',T1);      % enthalpy of N2 at T1 (kJ/kmol)
h_R = h_fuel+a*h_O2_R+(a*3.76)*h_N2_R;                 % total enthalpies of reactants (kJ/kmol)

% Products
h_CO2 = enthalpycalculator2('CO2',T2);      % enthalpy of CO2 at T2 (kJ/kmol)
h_H2O = enthalpycalculator2('H2O',T2);      % enthalpy of H2O at T2 (kJ/kmol)
h_O2_P = enthalpycalculator2('O2',T2);      % enthalpy of O2 at T2 (kJ/kmol)
h_N2_P = enthalpycalculator2('N2',T2);      % enthalpy of N2 at T2 (kJ/kmol)
h_P = x*h_CO2+(y/2)*h_H2O+(a-x-y/4)*h_O2_P+(a*3.76)*h_N2_P;            % total enthalpies of products (kJ/kmol)

% Moledular Weights
MW_C = 12;  % Molecular weight of carbon
MW_O = 16;  % Molecular weight of oxygen
MW_H = 1;   % Molecular weight of hydrogen
MW_N = 14;  % Molecular weight of nitrogen
n_MW = x*(MW_C+MW_O*2)+(y/2)*(MW_H*2+MW_O)+(a-x-(y/4))*(MW_O*2)+(a*3.76)*(MW_N*2);     % sum of number of moles of each product species times its molecular weight

delta_hc_kg = abs(h_R-h_P)/MW_fuel;        % lower heating value (kJ/kmol of fuel)
delta_hc = delta_hc_kg/((1/MW_fuel)*n_MW);

% Q3: Ratio of Pressure Increase and Volume Increase if the Max Pressure is
% P = 5.5 MPa

% Max Pressure
P3 = 5500; %kPa

%Ratio of Pressure Increase, alpha
alpha = P3/P2;

%Ratio of Volume Increase, beta

%Calculating Cp and Cv

%Givens/Known Values

T3 = alpha*T2;
R = 0.2891; %kPa*m^3/kg*K
V1 = (R*T1)/P1; %m^3/kg
V2 = V1/CR; %m^3/kg
V5 = V1;
V3 = V2;
P4 = P3;


%Cpbar & Cvbar (Calculated using specificheats.m)
cpbaro2 = 1.1491E4;
cvbaro2 = 1.1482E4;
cpbarN2 = 1.1798E4;
cvbarN2 = 1.1789E4;
cpbarfuel = 1.6374E5; %Adjust for fuel, recalculate using specificheats.m
cvbarfuel = 1.6373E5; %Adjust for fuel, recalculate using specificheats.m

cpmix = (cpbarfuel+(a*cpbaro2)+(3.76*a*cpbarN2))/(n_MW);
cvmix = (cvbarfuel+(a*cvbaro2)+(3.76*a*cvbarN2))/(n_MW);

beta = 1 + (delta_hc - cvmix*(alpha-1)*T2)/(cpmix*alpha*T2);


% Net work of the cycle
ExR = CR/beta;      % expansion ratio V5/V4
W_cycle = (P1*V1*CR^(n1-1))*((alpha*(beta-1)+((alpha*beta)/(n2-1))*(1-(1/ExR^(n2-1)))-(1/(n1-1))*(1-(1/(CR^n1-1)))));

% mep of cycle
mep = W_cycle/(V1-V2);

% thermal efficiency of cycle
eta_th = W_cycle/delta_hc;

% Kg of fuel per sec to produce 55 kW of power
Output= 55;
Input=(eta_th*Output)/100;

Fuel_rate=(Input/delta_hc);



function h=enthalpycalculatorfinal(Species,T)
if strcmpi(Species,'C8.26H15.5')==1
    a1=4.453623;
    a2=.03140168;
    a3=-0.12784105E-5;
    a4=0.02393996E-8;
    a5=-0.16690333E-13;
    a6=-0.04896696E6;

    theta=T/1000;
    h=4184*(a1*theta+a2*theta^2/2+a3*theta^3/3+a4*theta^4/4-a5*theta^-1+a6);
elseif strcmpi(Species,'C7.76H13.1')
    a1=-22.501;
    a2=227.99;
    a3=-177.26;
    a4=56.048;
    a5=0.4845;
    a6=-17.578;
    
    theta=T/1000;
    h=4184*(a1*theta+a2*theta^2/2+a3*theta^3/3+a4*theta^4/4-a5*theta^-1+a6);
elseif strcmpi(Species,'hydrogen')==1
    if T>=1000
        a1=0.02500000E2;
        a2=0;
        a3=0;
        a4=0;
        a5=0;
        a6=0.02547162E6;
        
        T_h=T;
        h=(a6+(a1*T_h)+((a2/2)*(T_h)^2)+((a3/3)*(T_h)^3)+((a4/4)*(T_h)^4)+((a5/5)*(T_h)^5))*8.314;
    end
    
end

end 
    
function h=enthalpycalculator2(Species,T)

if strcmpi(Species,'CO2')==1
    if T>=1000
        a1=0.04453623E2;
        a2=0.03140168E-1;
        a3=-0.12784105E-5;
        a4=0.02393996E-8;
        a5=-0.16690333E-13;
        a6=-0.04896696E6;
    elseif T<1000 
        a1=0.02275724E2;
        a2=0.09922072E-1;
        a3=-0.10409113E-4;
        a4=0.06866686E-7;
        a5=-0.02117280E-10;
        a6=-0.04837314E6;
    end
        
elseif strcmpi(Species,'O2')==1
    if T>=1000
        a1=0.03697578E2;
        a2=0.06135197E-2;
        a3=-0.12588420E-6;
        a4=0.01775281E-9;
        a5=-0.11364354E-14;
        a6=-0.12339301E4;
    elseif T<1000
        a1=0.03212936E2;
        a2=0.11274864E-2;
        a3=-0.05756150E-5;
        a4=0.13138773E-8;
        a5=-0.08768554E-11;
        a6=-0.10052490E4;
    end

        
elseif strcmpi(Species,'N2')==1
    if T>=1000
        a1=0.02926640E2;
        a2=0.14879768E-2;
        a3=-0.05684760E-5;
        a4=0.10097038E-9;
        a5=-0.06753351E-13;
        a6=-0.09227977E4;
    elseif T<1000
        a1=0.03298677E2;
        a2=0.14082404E-2;
        a3=-0.03963222E-4;
        a4=0.05641515E-7;
        a5=-0.02444854E-10;
        a6=-0.10208999E4;
    end
        
elseif strcmpi(Species,'H2O')==1
    if T>=1000
        a1=0.02672145E02;
        a2=0.03056293E-1;
        a3=-0.08730260E-5;
        a4=0.12009964E-9;
        a5=-0.06391618E-13;
        a6=-0.02989921E6;
    elseif T<1000
        a1=0.03386842E2;
        a2=0.03474982E-1;
        a3=-0.06354696E-4;
        a4=0.06968581E-7;
        a5=-0.02506588E-10;
        a6=-0.03020811E6;
        
    end
end


     h=8.314*(a6+a1*T+(a2/2)*T^2+(a3/3)*T^3+(a4/4)*T^4+(a5/5)*T^5);
     
end 

function [a,b,c,d,x,y,MW_fuel]= atombalancefinal(Fuel, ER) 

if strcmpi(Fuel,'C8.26H15.5')==1 
    MW_fuel=114.8;
    x=8.26;
    y=15.5; 
    
elseif strcmpi(Fuel,'C7.76H13.1')==1 
    MW_fuel=106.4;
    x=7.76;
    y=13.1; 
    
elseif strcmpi(Fuel,'H2')==1 
    MW_fuel=2;
    x=0;
    y=1;
    
end

astoich=x+ y/4;
a=astoich./ER;
b=x;
c=y/2;
d=a-b-(c/2);

end

%just type specificheats in the command window and it should start running
function [cpbar, cvbar]=specificheats

prompt1='What is your temperature(K): ';
T=input(prompt1);
theta=T/1000;
prompt2='What is your Fuel(Molecular Formula): ';
Fuel=input(prompt2,'s');
if strcmpi(Fuel,'C8.26H15.5')==1
    Mw=114.8;
    a1=-24.078;
    a2=256.63;
    a3=-201.68;
    a4=64.750;
    a5=0.5808;
elseif strcmpi(Fuel,'C7.76H13.1')
    Mw=106.4;
    a1=-22.501;
    a2=227.99;
    a3=-177.26;
    a4=56.048;
    a5=0.4845;
elseif T>1000 && strcmpi(Fuel,'H2')
    Mw=2;
    a1=2.991423;
    a2=0.0007000644;
    a3=-0.05633828E-6;
    a4=-0.09231578E-10;
    a5=0.15827519E-14;
elseif T<1000 && strcmpi(Fuel,'H2')
    Mw=2;
    a1=3.298124;
    a2=0.0008249441;
    a3=-0.08143015E-5;
    a4=-0.09475434E-9;
    a5=0.04134872E-11;
elseif T<1000 && strcmpi(Fuel, 'O2')
    Mw = 31.998;
    a1 = 0.03212936E2;
    a2 = 0.1127486E-2;
    a3 = -0.05756150E-5;
    a4 = 0.13138773E-8;
    a5 = -0.08768554E-11;
elseif T<1000 && strcmpi(Fuel, 'N2')
    Mw = 28.014;
    a1 = 0.03298677E2;
    a2 = 0.14082404E-2;
    a3 = -0.03963222E-4;
    a4 = 0.05641515E-7;
    a5 = -0.02444854E-10;
end


cpbar=4184*((a1*theta)+((a2*theta^2)/2)+((a3*theta^3)/3)+((a4*theta^4)/4)-(a5*theta^-1))
cp=cpbar/Mw;
R=8.314/Mw;
cv=cp-R;
cvbar=cv*Mw



end