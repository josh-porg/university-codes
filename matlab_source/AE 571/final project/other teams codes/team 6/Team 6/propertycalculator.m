function [MW,h,cp] = propertycalculator(Species,Temp)

% Constants
Ru = 8.3145;

% Reactants ---------------------------------------------------------------
% Methane
if strcmpi(Species,'CH4') == 1
    a1 = -0.29149;
    a2 = 26.327;
    a3 = -10.610;
    a4 = 1.5656;
    a5 = 0.16753;
    a6 = -18.331;
% Reactants cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 16.043;
h = 4184*(a1*(T/1000) + a2*((T/1000)^2)/2 + a3*((T/1000)^3)/3 + a4*((T/1000)^4)/4 - a5*((T/1000)^-1) + a6);


% Propane
elseif strcmpi(Species,'C3H8') == 1
    a1 = -1.4867;
    a2 = 74.339;
    a3 = -36.065;
    a4 = 8.0543;
    a5 = 0.01219;
    a6 = -27.313;
% Reactants cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 44.096;
h = 4184*(a1*(T/1000) + a2*((T/1000)^2)/2 + a3*((T/1000)^3)/3 + a4*((T/1000)^4)/4 - a5*((T/1000)^-1) + a6);

% Hexane
elseif strcmpi(Species,'C6H14') == 1
    a1 = -20.777;
    a2 = 210.48;
    a3 = -164.125;
    a4 = 52.832;
    a5 = 0.56635;
    a6 = -39.836;
% Reactants cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 86.177;
h = 4184*(a1*(T/1000) + a2*((T/1000)^2)/2 + a3*((T/1000)^3)/3 + a4*((T/1000)^4)/4 - a5*((T/1000)^-1) + a6);

% Isooctane
elseif strcmpi(Species,'C8H18') == 1
    a1 = -0.55313;
    a2 = 181.62;
    a3 = -97.787;
    a4 = 20.402;
    a5 = -0.03095;
    a6 = -60.751;
% Reactants cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 114.23;
h = 4184*(a1*(Tf/1000) + a2*((T/1000)^2)/2 + a3*((T/1000)^3)/3 + a4*((T/1000)^4)/4 - a5*((T/1000)^-1) + a6);

% Gasoline
elseif strcmpi(Species,'C8.26H15.5') == 1
    a1 = -24.078;
    a2 = 256.63;
    a3 = -201.68;
    a4 = 64.750;
    a5 = 0.5808;
    a6 = -27.562;
% Reactants cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 114.62;
h = 4184*(a1*(T/1000) + a2*((T/1000)^2)/2 + a3*((T/1000)^3)/3 + a4*((T/1000)^4)/4 - a5*((T/1000)^-1) + a6);

% Diesel
elseif strcmpi(Species,'C10.8H18.7') == 1
    a1 = -9.1063;
    a2 = 246.97;
    a3 = -143.74;
    a4 = 32.329;
    a5 = 0.0518;
    a6 = -50.128;
% Reactants cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 148.6;
h = 4184*(a1*(T/1000) + a2*((T/1000)^2)/2 + a3*((T/1000)^3)/3 + a4*((T/1000)^4)/4 - a5*((T/1000)^-1) + a6);

% Products ----------------------------------------------------------------
elseif strcmpi(Species,'CO') == 1
    if Temp >= 1000
        a1 = 0.03025078E2;
        a2 = 0.14426885E-2;
        a3 = -0.05630827E-5;
        a4 = 0.10185813E-9;
        a5 = -0.06910951E-13;
        a6 = -0.14268350E5;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 28;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    elseif Temp < 1000
        a1 = 0.03262451E2;
        a2 = 0.15119409E-2;
        a3 = -0.03881755E-4;
        a4 = 0.05581944E-7;
        a5 = -0.02474951E-10;
        a6 = -0.14310539E5;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 28;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    end
elseif strcmpi(Species,'CO2') == 1
    if Temp >= 1000
        a1 = 0.04453623E2;
        a2 = 0.03140168E-1;
        a3 = -0.12784105E-5;
        a4 = 0.02393996E-8;
        a5 = -0.16690333E-13;
        a6 = -0.04896696E6;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 44;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    elseif Temp < 1000
        a1 = 0.022725724E2;
        a2 = 0.09922072E-1;
        a3 = -0.10409113E-4;
        a4 = 0.06866686E-7;
        a5 = -0.02117280E-10;
        a6 = -0.04837314E6;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 44;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    end
elseif strcmpi(Species,'H2O') == 1
    if Temp >= 1000
        a1 = 0.02672145E2;
        a2 = 0.03056293E-1;
        a3 = -0.08730260E-5;
        a4 = 0.12009964E-9;
        a5 = -0.06391618E-13;
        a6 = -0.02989921E6;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 18;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    elseif Temp < 1000
        a1 = 0.03386842E2;
        a2 = 0.03474982E-1;
        a3 = -0.06354696E-4;
        a4 = 0.06968581E-7;
        a5 = -0.02506588E-10;
        a6 = -0.03020811E6;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 18;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    end
elseif strcmpi(Species,'O2') == 1
    if Temp >= 1000
        a1 = 0.03697578E2;
        a2 = 0.06135197E-2;
        a3 = -0.12588420E-6;
        a4 = 0.01775281E-9;
        a5 = -0.11364354E-14;
        a6 = -0.12339301E4;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 32;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    elseif Temp < 1000
        a1 = 0.03212936E2;
        a2 = 0.11274864E-2;
        a3 = -0.05756150E-5;
        a4 = 0.13138773E-8;
        a5 = -0.08768554E-11;
        a6 = -0.10052490E4;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 32;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    end
elseif strcmpi(Species,'N2') == 1
    if Temp >= 1000
        a1 = 0.02926640E2;
        a2 = 0.14879768E-2;
        a3 = -0.05684760E-5;
        a4 = 0.10097038E-9;
        a5 = -0.06753351E-13;
        a6 = -0.09227977E4;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 28;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    elseif Temp < 1000
        a1 = 0.03298677E2;
        a2 = 0.14082404E-2;
        a3 = -0.03963222E-4;
        a4 = 0.05641515E-7;
        a5 = -0.02444854E-10;
        a6 = -0.10208999E4;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 28;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    end
elseif strcmpi(Species,'H2') == 1
    if Temp >= 1000
        a1 = 0.02991423E2;
        a2 = 0.07000644E-2;
        a3 = -0.05633828E-6;
        a4 = -0.09231578E-10;
        a5 = 0.15827519E-14;
        a6 = -0.08350340E4;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 2;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    elseif Temp < 1000
        a1 = 0.03298124E2;
        a2 = 0.08249441E-2;
        a3 = -0.08143015E-5;
        a4 = -0.09475434E-9;
        a5 = 0.04134872E-11;
        a6 = -0.10125209E4;
% Product cp Calculator
syms T
theta = T/1000;
cp = 4.184*(a1 + a2*(theta) + a3*(theta^2) + a4*(theta^3) + a5*(theta^-2));
MW = 2;
h = Ru*(a6 + a1*T + (a2/2)*T^2 + (a3/3)*T^3 + (a4/4)*T^4 + (a5/5)*T^5);

    end

end