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