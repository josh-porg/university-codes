%AE571 Part I
clc;
clear;

C_Ratio = 16;
E_Ratio = 0.8;
T = 298;
P1 = 100;
n1 = 1.38;
n2 = 1.25;
Ru = 8.3145;

%Del1 = C_Ratio/beta1;

%hbar = Ru*(a6+a1*T1+(a2/2)*T1^2+(a3/3)*T^3+(a4/4)*T1^4+(a5/5)*T1^5);

% Part 1

% Dual Cycle [if gasoline is used as fuel]

% Four-Stroke Compression Ignition

% Q1) Find-Chemical Reaction Equation


% Gasoline = C8.26 H15.5

MW1 = 114.8;

x1 = 8.26;
y1 = 15.5;
z1 = (x1+y1/4)/E_Ratio; % usually a in the professors work

Chemical_Equation = 'C8.26H15.5 + 15.1687*(O2 + 3.76N2) --> 8.26*CO2 + 7.5*H2O + 3.0337*O2 + 57.0345*N2';

% number of moles for products

NO = 8.26*2+(15.5/2)+(z1-8.26-15.5/4)*2;
NC = 8.26;
NH = (15.5/2);
NN = 3.76*z1*2;


Species = "Gasoline";
% Q2) Lower heat value of the fuel at T2 (the temperature at the end of compression)

% Enthalpy inputs

a1=-24.078;                     %fuel enthalpy inputs
a2=256.63;                      %fuel enthalpy inputs
a3=-201.68;                     %fuel enthalpy inputs
a4=64.75;                       %fuel enthalpy inputs
a5=0.5808;                      %fuel enthalpy inputs
a6=-27.562;                     %fuel enthalpy inputs

O1 = 0.03212936*10^2;           % O2 temp for 300-1000  enthalpy inputs
O2 = 0.11274864*10^-2;          % O2 temp for 300-1000  enthalpy inputs
O3 = -0.0576150*10^-5;          % O2 temp for 300-1000  enthalpy inputs
O4 = 0.1318773*10^-8;           % O2 temp for 300-1000  enthalpy inputs
O5 = -0.0876855*10^-11;         % O2 temp for 300-1000  enthalpy inputs
O6 = -0.10052490*10^4;          % O2 temp for 300-1000  enthalpy inputs

N1 = 0.03298677*10^2;           % N2 temp 300-1000  enthalpy inputs
N2 = 0.14082404*10^-2;          % N2 temp 300-1000  enthalpy inputs
N3 = -0.03963222*10^-4;         % N2 temp 300-1000  enthalpy inputs
N4 = 0.05641515*10^-7;          % N2 temp 300-1000  enthalpy inputs
N5 = -0.02444854*10^-10;        % N2 temp 300-1000  enthalpy inputs
N6 = -0.10208999*10^4;          % N2 temp 300-1000  enthalpy inputs

CO2_1 = 0.02275724*(10^2);      % CO2 a1 for 300-1000  enthalpy inputs
CO2_2 = 0.09922072*(10^-1);     % CO2 a2 for 300-1000  enthalpy inputs
CO2_3 = -0.10409113*(10^-4);    % CO2 a3 for 300-1000  enthalpy inputs
CO2_4 = 0.06866686*(10^-7);     % CO2 a4 for 300-1000  enthalpy inputs
CO2_5 = -0.02117280*(10^-10);   % CO2 a5 for 300-1000  enthalpy inputs
CO2_6 = -0.04837314*(10^6);     % CO2 a6 for 300-1000  enthalpy inputs
  
H2O_1 = 0.03386842*(10^2);      % H2O a1 for 300-1000  enthalpy inputs
H2O_2 = 0.03474982*(10^-1);     % H2O a2 for 300-1000  enthalpy inputs
H2O_3 = -0.06354696*(10^-4);    % H2O a3 for 300-1000  enthalpy inputs
H2O_4 = 0.06968581*(10^-7);     % H2O a4 for 300-1000  enthalpy inputs
H2O_5 = -0.02506588*(10^-10);   % H2O a5 for 300-1000  enthalpy inputs
H2O_6 = -0.03020811*(10^6);     % H2O a6 for 300-1000  enthalpy inputs

O2_1 = 0.03212936*(10^2);       % O2 a1 for 300-1000  enthalpy inputs
O2_2 = 0.11274864*(10^-2);      % O2 a2 for 300-1000  enthalpy inputs
O2_3 = -0.05756150*(10^-5);     % O2 a3 for 300-1000  enthalpy inputs
O2_4 = 0.13138773*(10^-8);      % O2 a4 for 300-1000  enthalpy inputs
O2_5 = -0.08768554*(10^-11);    % O2 a5 for 300-1000  enthalpy inputs
O2_6 = -0.10052490*(10^4);      % O2 a6 for 300-1000  enthalpy inputs

N2_1 = 0.03298677*(10^2);       % N2 a1 for 300-1000  enthalpy inputs
N2_2 = 0.14082404*(10^-2);      % N2 a2 for 300-1000  enthalpy inputs
N2_3 = -0.03963222*(10^-4);     % N2 a3 for 300-1000  enthalpy inputs
N2_4 = 0.05641515*(10^-7);      % N2 a4 for 300-1000  enthalpy inputs
N2_5 = -0.02444854*(10^-10);    % N2 a5 for 300-1000  enthalpy inputs
N2_6 = -0.10208999*(10^4);      % N2 a6 for 300-1000  enthalpy inputs

O2T3_1 = 3.697578;
O2T3_2 = 0.0006135197;
O2T3_3 = -0.00000012588420;
O2T3_4 = 0.00000000001775281;
O2T3_5 = -0.0000000000000011364354;
O2T3_6 = -1233.9301;

N2T3_1 = 2.926640;
N2T3_2 = 0.0014879768;
N2T3_3 = -0.0000005684760;
N2T3_4 = 0.00000000010097038;
N2T3_5 = -0.000000000000006753351;
N2T3_6 = -922.7977;

CO2T3_1 = 4.453623;
CO2T3_2 = 0.003140168;
CO2T3_3 = -0.0000012784105;
CO2T3_4 = 0.0000000002393996;
CO2T3_5 = -0.000000000000016690333;
CO2T3_6 = -48966.96;

H2OT3_1 = 2.672145;
H2OT3_2 = 0.003056293;
H2OT3_3 = -0.0000008730260;
H2OT3_4 = 0.00000000012009964;
H2OT3_5 = -0.000000000000006391618;
H2OT3_6 = -29899.21;
   
% molecular weights

MWO = 15.999;                   % Molecular Weight of Oxygen
MWC = 12.011;                    % Molecular Weight of Carbon
MWN = 14.007;                   % Molecular Weight of Nitrogen
MWH = 1.008;                    % Molecular Weight of Hydrogen

% T2 for LHV

% LHV

Theta = 298/1000;

T2 = T*C_Ratio^(n1-1);
Theta2 = T2/1000;

H1 = 4184*(a1*Theta2+((a2)/2)*((Theta2)^2)+((a3)/3)*((Theta2)^3)+((a4)/4)*((Theta2)^4)+a5/5*((Theta2)^5)+a6);            % fuel reactant enthalpy
H2 = Ru*(O1*T2+(O2/2)*(T2^2)+(O3/3)*((T2)^3)+(O4/4)*((T2)^4)+(O5/5)*((T2)^5)+(O6));                                 % O2 reactant enthalpy
H3 = Ru*(N1*T2+(N2/2)*(T2^2)+(N3/3)*((T2)^3)+(N4/4)*((T2)^4)+(N5/5)*((T2)^5)+(N6));                                 % N2 reactant enthalpy
H4 = Ru*(CO2_1*T2+((CO2_2)/2)*((T2)^2)+((CO2_3)/3)*((T2)^3)+((CO2_4)/4)*((T2)^4)+((CO2_5)/5)*((T2)^5)+(CO2_6));     % CO2 product enthalpy
H5 = Ru*(H2O_1)*T2+((H2O_2)/2)*((T2)^2)+((H2O_3)/3)*((T2)^3)+((H2O_4)/4)*((T2)^4)+((H2O_5)/5)*((T2)^5)+(H2O_6);     % H2O product enthalpy
H6 = Ru*(O2_1)*T2+((O2_2)/2)*((T2)^2)+((O2_3)/3)*((T2)^3)+((O2_4)/4)*((T2)^4)+((O2_5)/5)*((T2)^5)+(O2_6);           % O2 product enthalpy
H7 = Ru*(N2_1)*T2+((N2_2)/2)*((T2)^2)+((N2_3)/3)*((T2)^3)+((N2_4)/4)*((T2)^4)+((N2_5)/5)*((T2)^5)+(N2_6);           % N2 product enthalpy

h_reactants1= (H1 + z1*H2 + z1*3.76*H3);                            % Total enthalpy of reactants
        
h_products1= (x1*H4 + (y1/2)*H5 + (z1-x1-y1/4)*H6 + z1*3.76*H7);    % Total enthalpy of products

Change_Hc1 = h_reactants1 - h_products1;

Change_hcbar1 = Change_Hc1/1;

Change_hc = Change_hcbar1/MW1;

unitchange1 = (MWO*(8.26*2+15.5/2+(z1-8.26-15.5/4)*2)+MWN*3.76*z1*2+MWC*8.26+MWH*2*15.5/2)/MW1;

LHV1 = Change_hc/unitchange1;
 

% Q3) Ratio of pressure increase Alpha and volume increase Beta if the
% maximum pressure the engine can sustain is 5.5 Mpa
P3 = 5500000;
alpha1 = P3/(P1*1000*C_Ratio^(n1));
T3 = T2*alpha1;
Theta3 = T3/1000;

%Cp1 = 4.184*(a1+a2*Theta+a3*(Theta^2)+a4*(Theta^3)+a5*(Theta^-2));                                                      % fuel reactant specific heat
%Cp2 = 4.184*(O1+O2*Theta+O3*(Theta^2)+O4*(Theta^3)+O5*(Theta^-2));                                                      % O2 reactant specific heat
%Cp3 = 4.184*(N1+N2*Theta+N3*(Theta^2)+N4*(Theta^3)+N5*(Theta^-2));                                                      % N2 reactant specific heat
CpCO2 = Ru*(CO2_1+CO2_2*T2+CO2_3*(T2^2)+CO2_4*(T2^3)+CO2_5*(T2^4));                                       % CO2 product specific heat
CpH20 = Ru*(H2O_1+H2O_2*T2+H2O_3*(T2^2)+H2O_4*(T2^3)+H2O_5*(T2^4));                                       % H2O product specific heat
CpO2 = Ru*(O2_1+O2_2*T2+O2_3*(T2^2)+O2_4*(T2^3)+O2_5*(T2^4));                                            % O2 product specific heat
CpN2 = Ru*(N2_1+N2_2*T2+N2_3*(T2^2)+N2_4*(T2^3)+N2_5*(T2^4));

CpCO2t3 = Ru*(CO2T3_1+CO2T3_2*T3+CO2T3_3*(T3^2)+CO2T3_4*(T3^3)+CO2T3_5*(T3^4));                                       % CO2 product specific heat
CpH20t3 = Ru*(H2OT3_1+H2OT3_2*T3+H2OT3_3*(T3^2)+H2OT3_4*(T3^3)+H2OT3_5*(T3^4));                                       % H2O product specific heat
CpO2t3 = Ru*(O2T3_1+O2T3_2*T3+O2T3_3*(T3^2)+O2T3_4*(T3^3)+O2T3_5*(T3^4));                                            % O2 product specific heat
CpN2t3 = Ru*(N2T3_1+N2T3_2*T3+N2T3_3*(T3^2)+N2T3_4*(T3^3)+N2T3_5*(T3^4)); 
% N2 product specific heat

Tguess1 = 1000;
Tref = 298;

for i = 1:4000

    Tmp(i) = Tguess1 + 0.75*i;

    [h0_fuelTref,h0_O2Tref,h0_N2Tref,h0_CO2Tref,h0_H2OTref,h0_fuelT1,h0_O2T1,h0_N2T1,cp_fuelT1,h0_CO2T4,h0_H2OT4,h0_O2T4,h0_N2T4] = AE571_Part1_Calculators(Species,Tmp,Tref,i);
% Enthalpies of Reactants   T1 = temp of reactants
    h_fuel = h0_fuelTref;
    h_O2T1 = h0_O2Tref;
    h_N2T1 = h0_N2Tref;

% Enthalpies of Products
    h_CO2 = h0_CO2T4;         % - h0_CO2Tref + (h0_CO2T2 - h0_CO2Tref);
    h_H2O = h0_H2OT4;         %h0_H2OTref + (h0_H2OT2 - h0_H2OTref);
    h_O2T2 = h0_O2T4;         %h0_O2Tref + (h0_O2T2 - h0_O2Tref);
    h_N2T2 = h0_N2T4;

    h_reactants = h_fuel + z1*(h_O2T1 + 3.76*h_N2T1);
    h_products = 8.26*h_CO2 + (15.5/2)*h_H2O + (z1-8.26-15.5/4)*h_O2T2 + 3.76*z1*h_N2T2;

error(i) = abs((h_products - h_reactants)/h_reactants);
%error(i) = abs((h_productsairER8 - h_reactantsairER8)/h_reactantsairER8);


    if error(i) < 0.01
        break
    end

end

hold on
figure(1)
title('T vs Iteration')
plot(1:i,Tmp) %i is iteration #

figure(2)
title('Error vs Iteration')
plot(1:i,error) %plot error vs iteration #

T4 = Tmp(end);
CpCO2t4 = Ru*(CO2T3_1+CO2T3_2*T4+CO2T3_3*(T4^2)+CO2T3_4*(T4^3)+CO2T3_5*(T4^4));                                       % CO2 product specific heat
CpH20t4 = Ru*(H2OT3_1+H2OT3_2*T4+H2OT3_3*(T4^2)+H2OT3_4*(T4^3)+H2OT3_5*(T4^4));                                       % H2O product specific heat
CpO2t4 = Ru*(O2T3_1+O2T3_2*T4+O2T3_3*(T4^2)+O2T3_4*(T4^3)+O2T3_5*(T4^4));                                            % O2 product specific heat
CpN2t4 = Ru*(N2T3_1+N2T3_2*T4+N2T3_3*(T4^2)+N2T3_4*(T4^3)+N2T3_5*(T4^4)); 


Cpmixt2 = (8.26*CpCO2+3.76*z1*CpN2+(15.5/2)*CpH20+(z1-8.26-15.5/4)*CpO2)/(MWO*(8.26*2+15.5/2+(z1-8.26-15.5/4)*2)+MWN*3.76*z1*2+MWC*8.26+MWH*2*15.5/2);

Cvmixt2 = (8.26*(CpCO2-Ru)+3.76*z1*(CpN2-Ru)+(15.5/2)*(CpH20-Ru)+(z1-8.26-15.5/4)*(CpO2-Ru))/(MWO*(8.26*2+15.5/2+(z1-8.26-15.5/4)*2)+MWN*3.76*z1*2+MWC*8.26+MWH*2*15.5/2);

Cpmixt3 = (8.26*CpCO2t3+3.76*z1*CpN2t3+(15.5/2)*CpH20t3+(z1-8.26-15.5/4)*CpO2t3)/(MWO*(8.26*2+15.5/2+(z1-8.26-15.5/4)*2)+MWN*3.76*z1*2+MWC*8.26+MWH*2*15.5/2);

Cvmixt3 = (8.26*(CpCO2t3-Ru)+3.76*z1*(CpN2t3-Ru)+(15.5/2)*(CpH20t3-Ru)+(z1-8.26-15.5/4)*(CpO2t3-Ru))/(MWO*(8.26*2+15.5/2+(z1-8.26-15.5/4)*2)+MWN*3.76*z1*2+MWC*8.26+MWH*2*15.5/2);

Cpmixt4 = (8.26*CpCO2t4+3.76*z1*CpN2t4+(15.5/2)*CpH20t4+(z1-8.26-15.5/4)*CpO2t4)/(MWO*(8.26*2+15.5/2+(z1-8.26-15.5/4)*2)+MWN*3.76*z1*2+MWC*8.26+MWH*2*15.5/2);

Cvmixt4 = (8.26*(CpCO2t4-Ru)+3.76*z1*(CpN2t4-Ru)+(15.5/2)*(CpH20t4-Ru)+(z1-8.26-15.5/4)*(CpO2t4-Ru))/(MWO*(8.26*2+15.5/2+(z1-8.26-15.5/4)*2)+MWN*3.76*z1*2+MWC*8.26+MWH*2*15.5/2);

beta1 = (LHV1 - (Cvmixt3*T3 - Cvmixt2*T2) + Cpmixt3*T3)/(T3*Cpmixt4);
% beta1 = 1+ ((LHV1 - (Cvmixt3*T3-Cvmixt2*T2))/(Cvmixt4*T3));

%beta1 = 2109.75/T3;
%beta1 = 1+(LHV1-Cvmix*(alpha1-1)*T2)/(Cpmix*alpha1*T2);

%R1=Ru/Mw1;
%alpha1 = (h_reactants1/(Cv1*T1)+C_Ratio^(n1-1))/C_Ratio^(n1-1);

% Q4) Net work of the cycle

del1 = C_Ratio/beta1;

MWreacant = (MWC*8.26+MWH*15.5+z1*(MWO*2+MWN*2*3.76))/(1+z1*(1+3.76));

V1 = ((Ru/MWreacant)*T)/P1;
W1 = P1*V1*C_Ratio^(n1-1)*(alpha1*(beta1-1)+(alpha1*beta1/(n2-1))*(1-(1/(del1^(n2-1))))-(1/(n1-1))*(1-(1/(C_Ratio^(n1-1)))));

% Q5) MEP of the cycle

V2 = V1/C_Ratio;

MEP1 = W1/(V1-V2);

% Q6) Thermal Efficiency of the cycle

Thermal_effeciency1 = W1/LHV1;

% Q7) How many kilograms of fuel per sec needs to be injected to the engine
% to produce 55 kW of power?

work_ratio = 55/W1;
Fuel_KG = work_ratio;

qdot = 55/Thermal_effeciency1;
Fuel_KG_REQ = qdot/LHV1;

wr = 55/W1;

colNames = {'Variable Names','Values'};
rowNames = {'LHV','alpha','beta','Work','MEP','Thermal Efficiency','Fuel Required'};
Outputval =[LHV1,"kJ/kg";alpha1,"N/A";beta1,"N/A";W1,"kJ/kg";MEP1,"kPa";Thermal_effeciency1,"N/A";Fuel_KG_REQ,"kg/s"];
%Outputval =[LHV1,0;alpha1,0;beta1,0;W1,0;MEP1,0;Thermal_effeciency1,0;Fuel_KG_REQ,0];
answers = array2table(Outputval,'RowNames',rowNames,'VariableNames',colNames);

disp(Chemical_Equation)

answers