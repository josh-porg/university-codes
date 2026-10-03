%Richard Bramlette
%Graduate Research Assistant
%University of Kansas
%27 July 2010

%PROGRAM TO DESIGN A GLUHAREFF PRESSURE JET ENGINES

%This program will allow the user to design and size a Gluhareff Pressure
%Jet engine based on the parameters available such as engine thrust, tank
%pressure, etc...

%Design Input Specifications
ID0 = 0.0400;                       %Injector Inner Diameter (in.)
Pf = 74.700;                        %Injector Feed Pressure (psia)
Tf = 420.00;                        %Injector Feed Temp (deg F)
Tflame = 1670.0;                    %Third Stage Flame Temp (deg F)

N0 = 0.9097;                        %Injector Efficiency (per 1)
Ni = 0.6747;                        %Inlet Efficiency (per 1)

phi3 = 0.9429;                      %Third Stage Equivalency Ratio
phi2 = 0.8068;                      %Second Stage Equivalency Ratio
phi1 = 0.3807;                      %First Stage Equivalency Ratio
f_stoich = 15.500;                  %Stoichiometric Air/Fuel Ratio for Fuel

F15 = 639;                          %G8-2-15 Operating Frequency (Hz)
ID15 = 0.125;                       %G8-2-15 Injector Inner Diameter (in.)

%Operating Conditions
Pamb = 14.700;                      %Ambient Pressure (psia)
Tamb = 70.000;                      %Ambient Temp (deg F)
g = 32.178;                         %Acceleration Due To Gravity (ft/sec^2)
y_prop = 1.15;                      %Ratio of Specific Heats for Propane
y_air = 1.40;                       %Ratio of Specific Heats for Air
R_prop = 1130;                      %Propane Gas Constant (ft-lbf/slug-R)
R_air = 1716;                       %Air Gas Constant (ft-lbf/slug-R)

% Derived Specifications
%========================
A0 = pi*(ID0/2)^2;                  %Injector Exit Area (sq.in.)
Tf = Tf + 459.67;                   %Injector Temperature (Rankine)
Tamb = Tamb + 459.67;               %Ambient Temperature (Rankine)
rhoAMB = Pamb*144/(R_air*Tamb);     %Ambient Density (slug/ft^3)
Tflame = Tflame + 459.67;           %Third Stage Flame Temp (Rankine)

%=========================================================================
%                        Injection Nozzle Analysis
%=========================================================================
%Injector Mass Flow Rate (slug/sec)
yf = y_prop; Rf = R_prop;
mdot0_ideal = Pf*A0*sqrt((yf/(Rf*Tf))*(2/(yf+1))^((yf+1)/(yf-1)));
mdot0 = N0*mdot0_ideal;

%Injector Exit Temperature (deg R)
y0 = yf;
T0 = Tf*(Pamb/Pf)^((y0-1)/y0);

%Injector Sound Speed (ft/sec)
R0 = Rf;
a0 = sqrt(y0*R0*T0);

%Injector Exit Density (slug/ft^3)
R0 = Rf; P0 = Pf;
rho0 = P0/(R0*T0);

%Injector Exit Speed (ft/sec)
V0 = N0*a0;

%Injector Dynamic Pressure (psia)
Q0 = (1/144)*(1/2)*rho0*V0^2;

%Injector Kinetic Energy (slug-ft^2/sec^2)
KE0 = (1/2)*mdot0*V0^2;

%Injector (and thus engine) Operating Frequency (Hz)
Fx = F15*(ID15/ID0);

%Injector Length (Unnecessary; Ignore)
L0 = 0;

%Save All Injector Data
Injector =  [0; N0; Fx; mdot0;   rho0;   L0; A0; ID0; V0; y0; R0;   T0-459.67; a0; P0; Q0; 1; KE0];
InjectorG = [0; N0; Fx; mdot0*g; rho0*g; L0; A0; ID0; V0; y0; R0/g; T0;        a0; P0; Q0; 1; KE0*g];

%=========================================================================
%                        Air Intake System Analysis
%=========================================================================

% Sizing FIRST Stage Inlet
%==========================
f1 = phi1*f_stoich;             %Necessary fuel/air mixture
mdot1 = f1*mdot0;               %Total mass flow for fuel/air mixture (slug/sec)

N1 = Ni;                        %Efficiency from Injector to First Stage (per 1)
KE1 = N1*KE0;                   %First Stage Kinetic Energy (slug-ft^2/sec^2)
V1 = sqrt(2*KE1/mdot1);         %First Stage Inlet Speed (ft/sec)

rho1 = rhoAMB*(1-(1/f1)) + rho0*(1/f1);     %Rule of Mixtures to find rho
A1 = 144*mdot1/(rho1*V1);       %First Stage Area needed to achieve mdot1 (in.^2)
ID1 = 2*sqrt(A1/pi);            %First Stage Inner Diameter needed for Area (in.)

Q1 = (1/144)*(1/2)*rho1*V1^2;   %First Stage Dynamic Pressure (psia)
P1 = Pamb - Q1;                 %First Stage Static Pressure (psia) ***CHECK THIS***

R1 = R_air*(1-(1/f1)) + R_prop*(1/f1);      %Rule of Mixtures to find R
y1 = y_air*(1-(1/f1)) + y_prop*(1/f1);      %Rule of Mixtures to find gamma
Tex = P1*144/(rho1*R1);                     %Ideal Gas Law Temp (R) ***CHECK THIS***
T1 = Tamb*(1-(1/f1)) + Tex*(1/f1);          %First Stage Temperature (R) ***CHECK THIS***

%R1 = Rf; y1 = yf; T1 = T0;     %Assuming Fuel Properties into First Stage
%R1 = R_air; y1 = y_air; T1 = Tamb;

a1 = sqrt(y1*R1*T1);            %Sound Speed in First Stage (ft/sec)
L1 = a1/(2*Fx);                 %First Stage Inlet Length (ft)

% Sizing SECOND Stage Inlet
%==========================
f2 = phi2*f_stoich;             %Necessary fuel/air mixture
mdot2 = f2*mdot0;               %Total mass flow for fuel/air mixture (slug/sec)

N2 = Ni^2;                      %Efficiency from Injector to Second Stage (per 1)
KE2 = N2*KE0;                   %Second Stage Kinetic Energy (slug-ft^2/sec^2)
V2 = sqrt(2*KE2/mdot2);         %Second Stage Inlet Speed (ft/sec)

rho2 = rhoAMB*(1-(1/f2)) + rho0*(1/f2);     %Rule of Mixtures to find rho
A2 = 144*mdot2/(rho2*V2);       %Second Stage Area needed to achieve mdot2 (in.^2)
ID2 = 2*sqrt(A2/pi);            %First Stage Inner Diameter needed for Area (in.)

Q2 = (1/144)*(1/2)*rho2*V2^2;   %Second Stage Dynamic Pressure (psia)
P2 = Pamb - Q2;                 %Second Stage Static Pressure (psia) ***CHECK THIS***

R2 = R_air*(1-(1/f2)) + R_prop*(1/f2);      %Rule of Mixtures to find R
y2 = y_air*(1-(1/f2)) + y_prop*(1/f2);      %Rule of Mixtures to find gamma
T2 = P2*144/(rho2*R2);          %Second Stage Temperature (R) ***CHECK THIS***

a2 = sqrt(y2*R2*T2);            %Sound Speed in First Stage (ft/sec)
L2 = a2/(2*Fx);                 %Second Stage Inlet Length (ft)

% Sizing THIRD Stage Inlet
%==========================
f3 = phi3*f_stoich;             %Necessary fuel/air mixture for combustion
mdot3 = f3*mdot0;               %Total mass flow for fuel/air mixture (slug/sec)

N3 = Ni^3;                      %Efficiency from Injector to Third Stage (per 1)
KE3 = N3*KE0;                   %Third Stage Kinetic Energy (slug-ft^2/sec^2)

V3 = sqrt(2*KE3/mdot3);         %Third Stage Inlet Speed (ft/sec)

rho3 = rhoAMB*(1-(1/f3)) + rho0*(1/f3);     %Rule of Mixtures to find rho
A3 = 144*mdot3/(rho3*V3);       %Third Stage Area needed to achieve mdot3 (in.^2)
ID3 = 2*sqrt(A3/pi);            %First Stage Inner Diameter needed for Area (in.)

Q3 = (1/144)*(1/2)*rho3*V3^2;   %Third Stage Dynamic Pressure (psia)
P3 = Pamb - Q3;                 %Third Stage Static Pressure (psia) ***CHECK THIS***

R3 = R_air*(1-(1/f3)) + R_prop*(1/f3);      %Rule of Mixtures to find R
y3 = y_air*(1-(1/f3)) + y_prop*(1/f3);      %Rule of Mixtures to find gamma
T3 = Tflame;                    %Third Stage Temperature (R) **ASSUMED**

a3 = sqrt(y3*R3*T3);            %Sound Speed in First Stage (ft/sec)
L3 = a3/(4*Fx);                 %Third Stage Inlet Length (ft)

%Save All Inlet Data
Inlet1 =  [1; N1; Fx; mdot1;   rho1;   L1; A1; ID1; V1; y1; R1;   T1-459.67; a1; P1; Q1;     f1; KE1];
Inlet2 =  [2; N2; Fx; mdot2;   rho2;   L2; A2; ID2; V2; y2; R2;   T2-459.67; a2; P2; Q2;     f2; KE2];
Inlet3 =  [3; N3; Fx; mdot3;   rho3;   L3; A3; ID3; V3; y3; R3;   T3-459.67; a3; P3; Q3;     f3; KE3];

Inlet1G = [1; N1; Fx; mdot1*g; rho1*g; L1; A1; ID1; V1; y1; R1/g; T1;        a1; P1; Q1/144; f1; KE1*g];
Inlet2G = [2; N2; Fx; mdot2*g; rho2*g; L2; A2; ID2; V2; y2; R2/g; T2;        a2; P2; Q2/144; f2; KE2*g];
Inlet3G = [3; N3; Fx; mdot3*g; rho3*g; L3; A3; ID3; V3; y3; R3/g; T3;        a3; P3; Q3/144; f3; KE3*g];

%=========================================================================
%                       Combustion Chamber Analysis
%=========================================================================

%Save All Combustor Data
%Combustor = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

%=========================================================================
%                           Exit Nozzle Analysis
%=========================================================================

%Save All Nozzle Data
%Nozzle = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

%=========================================================================
%                            Trend-Based Sizing
%=========================================================================

Ln = L2*2;                              %Length of Nozzle/L2 = 2 (ft)
IDcc = L3/2;                            %Diameter of Combustor/L3 = 1/2 (in.)
Lcc = ((3.6988*ID0) + 0.8307)*L2*1.15;  %Length of Combustor (ft)
            %Note: Trend of ratio with injector ID and adding 15% for nose
            %section (which was not included in trend calc).
IDn = (2/3)*IDcc;                       %IDn/IDcc = 0.66 (in.)
ID3f = ID3*(1.25);                      %Flared slope of 8 (in.)

%=========================================================================
%                             Organize Results
%=========================================================================

%Concatenate All Data
DATA = [Injector Inlet1 Inlet2 Inlet3];
DATA_G = [InjectorG Inlet1G Inlet2G Inlet3G];

%Write out designed specifications
fprintf('The Engine has been sized based on:\n')
fprintf('Injector Inner Diameter = %2.4f in.\n',ID0)
fprintf('Injector Feed Press.    = %4.1f psi\n',Pf)
fprintf('Injector Feed Temp.     = %4.1f deg F\n',Tf-459.67)
fprintf('Injector Efficiency     = %3.2f %%\n',N0*100)
fprintf('Inlet Efficiency        = %3.2f %%\n',Ni*100)
fprintf('Equivalency Ratios      = %1.4f %1.4f %1.4f\n',phi1,phi2,phi3)
fprintf('Stoich. Air/Fuel Ratio  = %3.2f\n',f_stoich)
fprintf('=============================================\n\n')
fprintf('Operating Frequency = %5.1f Hz\n\n',Fx)
fprintf('First Stage Inlet:\n')
fprintf('===================\n')
fprintf('Inner Diameter = %2.4f in.\n',ID1)
fprintf('Length         = %2.4f in.\n\n',L1*12)
fprintf('Second Stage Inlet:\n')
fprintf('===================\n')
fprintf('Inner Diameter = %2.4f in.\n',ID2)
fprintf('Length         = %2.4f in.\n\n',L2*12)
fprintf('Third Stage Inlet:\n')
fprintf('===================\n')
fprintf('Inner Diameter = %2.4f in.\n',ID3)
fprintf('Flared ID      = %2.4f in.\n',ID3f)
fprintf('Length         = %2.4f in.\n\n',L3*12)
fprintf('Combustion Chamber:\n')
fprintf('===================\n')
fprintf('Inner Diameter = %2.4f in.\n',IDcc*12)
fprintf('Length (w nose)= %2.4f in.\n\n',Lcc*12)
fprintf('Nozzle:\n')
fprintf('===================\n')
fprintf('Inner Diameter = %2.4f in.\n',IDn*12)
fprintf('Length         = %2.4f in.\n\n',Ln*12)
