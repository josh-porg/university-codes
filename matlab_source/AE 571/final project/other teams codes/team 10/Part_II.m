clc; 
clear;
%% Part 1
x = 0;
y = 2;
equivalence = .8;
astoich = x + y/4;
a = astoich/equivalence;
b = x;
c = y/2;
d = a - b - c/2;
an = 3.76*a;

%% Part 2
T1 = 298;
E = 16;
P1 = 100000;
n1 = 1.38;
n2 = 1.25;
T2 = T1*E^(n1-1);
theta = T2/1000;
hH2 =  8.3145*(-.10125209*10^4 + .03298124*10^2*T2 + .08249441*10^-2*.5*T2^2 - (1/3)*.08143015*10^-5*T2^3 - (1/4)*.09475434*10^-9*T2^4 + (1/5)*.04134872*10^-11*T2^5);
hO2 = 8.3145*(-.1005249E4 + .03212936E2*T2 + .11274864E-2*.5*T2^2 - (1/3)*.0575615E-5*T2^3 + (1/4)*.13138773E-8*T2^4 - (1/5)*.08768554E-11*T2^5);             
hN2 = 8.3145*(-.10208999E4 + .03298677E2*T2 + .14082404E-2*.5*T2^2 - (1/3)*.03963222E-4*T2^3 + (1/4)*.05641515E-7*T2^4 - (1/5)*.02444854E-10*T2^5); 
hH2O = 8.3145*(-.03020811E6 +.03386841E2*T2 + .03474982E-1*.5*T2^2 - (1/3)*.06354696E-4*T2^3 + (1/4)*.06968581E-7*T2^4 - (1/5)*.02506588E-10*T2^5);
HR = hH2 + a*hO2 + a*3.76*hN2;
HP = c*hH2O + d*hO2 + an*hN2;
HC = HR - HP;
MWH2 = 1.00784*2;
hcf = HC/MWH2;
mix = (c*18.0153 + d*32 + an*28.014)/MWH2;
hc = hcf/mix;

%% Part 3
P2 = P1*E^n1;
prompt = "What is the maximum sustainable pressure in MPa?\n";
mP3 = input(prompt);
P3 = mP3*1000*1000;
alpha = P3/P2;
T3 = T2*alpha;
beta = 1;
i = 1;
heat = 1;
while heat < hc
    T4 = beta*T3;
    Ru = 8.3145;
    cpH2 = @(X) Ru.*(.02991423E2 + .07000644E-2.*X + -.05633828E-6.*X.^2 - .09231578E-10.*X.^3 + .15827519E-14.*X.^4);
    cpO2 = @(Z) Ru.*(.03697578*10^2 + .06135197*10^-2.*Z - .12588420*10^-6.*Z.^2 + .01775281*10^-9.*Z.^3 + -.11364354*10^-14.*Z.^4);
    cpN2 = @(X) Ru.*(.0292664*10^2 + .14879768*10^-2.*X + -.0568476*10^-5.*X.^2 + .10097038*10^-9.*X.^3 + -.06753351*10^-13.*X.^4);
    CpH2 = integral(cpH2,T3,T4);
    CpO2 = integral(cpO2,T3,T4);
    CpN2 = integral(cpN2,T3,T4);
    cvO21 = @(X) Ru.*(.03212936*10^2 + .11274864*10^-2.*X - .05756150*10^-5.*X.^2 + .1318773*10^-8.*X.^3 + -.08768554*10^-11.*X.^4);
    cvN21 = @(X) Ru.*(.03298677E2 + .14082404E-2.*X + -.03963222E-4.*X.^2 + .05641515E-7.*X.^3 + -.02444854E-10.*X.^4);
    cvO22 = @(Z) Ru.*(.03697578*10^2 + .06135197*10^-2.*Z - .12588420*10^-6.*Z.^2 + .01775281*10^-9.*Z.^3 + -.11364354*10^-14.*Z.^4);
    cvN22 = @(X) Ru.*(.0292664*10^2 + .14879768*10^-2.*X + -.0568476*10^-5.*X.^2 + .10097038*10^-9.*X.^3 + -.06753351*10^-13.*X.^4);
    cvH21 = @(X) Ru.*(.03298124*10^2 + .08249441*10^-2.*X + (-.08143015*10^-5).*X.^2 - .09475434*10^-9.*X.^3 + .04134872*10^-11.*X.^4);
    cvH22 = @(X) Ru.*(.02991423E2 + .07000644E-2.*X + -.05633828E-6.*X.^2 - .09231578E-10.*X.^3 + .15827519E-14.*X.^4);
    CvO21 = integral(cvO21,T2,1000);
    CvN21 = integral(cvN21,T2,1000);
    CvO22 = integral(cvO22,1000,T3);
    CvN22 = integral(cvN22,1000,T3);
    CvH21 = integral(cvH21,T2,1000);
    CvH22 = integral(cvH22,1000,T3);
    CvH2 = CvH21 + CvH22;
    CvO2 = CvO21 + CvO22;
    CvN2 = CvN21 + CvN22;
    Cp = (CpH2 + a*CpO2 + an*CpN2)/(MWH2 + a*32 + an*28);
    Cv = (CvH2-Ru + a*(CvO2-Ru) + an*(CvN2-Ru))/(MWH2 + a*32 + an*28);
    heat = Cv + Cp;
    beta = beta + .001;
end

%% Part 4
delta = E/beta;
MWair = 28.9647;
R1 = Ru/MWair;
V1 = R1*1000*T1/P1;
wcycle = (P1*V1*E^(n1-1))*[alpha*(beta-1)+[alpha*beta/(n2-1)]*(1-1/delta^(n2-1))-(1/(n1-1))*(1-1/E^(n1-1))];
Wcycle = wcycle/1000;

%% Part 5
V2 = V1/E;
mep = wcycle/(V1-V2);
Mep = mep/1000;

%% Part 6
eta = Wcycle/hc;

%% Part 7
kgs = 55/hcf;



%% Print
clc;
formatSpec = '1) H2 + %3.3f(O2 + 3.76N2) = %1.1fH2O + %1.3fO2 + %5.2fN2.\n';
fprintf(formatSpec,a,c,d,an)
formatSpec = '2) The lower heating value for T2 is %6.2f KJ/kg.\n';
fprintf(formatSpec, hc)
formatSpec = '3) The ratio of pressure increase is %5.3f and the volume increase ratio is %4.3f.\n';
fprintf(formatSpec,alpha,beta)
formatSpec = '4) The net cycle work is %6.2f KJ/kg.\n';
fprintf(formatSpec,Wcycle)
formatSpec = '5) The mean effective pressure is %6.2f KPa.\n';
fprintf(formatSpec,Mep)
formatSpec = '6) The thermal efficiency of the cycle is %5.4f.\n';
fprintf(formatSpec,eta)
formatSpec = '7) %4.6f kg of fuel/sec is required to generate 55kW of power.\n';
fprintf(formatSpec,kgs)


