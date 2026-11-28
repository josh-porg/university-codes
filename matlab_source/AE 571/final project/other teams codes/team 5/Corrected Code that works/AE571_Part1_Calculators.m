function [h0_fuelTref,h0_O2Tref,h0_N2Tref,h0_CO2Tref,h0_H2OTref,h0_fuelT1,h0_O2T1,h0_N2T1,cp_fuelT1,h0_CO2T4,h0_H2OT4,h0_O2T4,h0_N2T4] = Copy_of_enthalpies(Species,Tmp,Tref,i)

if strcmpi(Species,'Gasoline') == 1
    %mw = 16.043;
    a1 = -24.078;
    a2 = 256.63;
    a3 = -201.68;
    a4 = 64.75;
    a5 = 0.5808;
    a6 = -27.562; 

end
%theta = Tmp/1000;
T = Tref;
thetaT1 = T/1000;
thetaTref = Tref/1000;

% Note that in this problem, Tref and T1 are equal to each other and are < 1000K, so they can both use the same coefficients
% Code will need to be altered if T1 is not equal to Tref, and also make a case for T1 > 1000K
% Also note that T2 = T

% O2
a1O2Tref = 3.212936;
a2O2Tref = 0.001127486;
a3O2Tref = -0.0000005756150;
a4O2Tref = 0.0000000013138773;
a5O2Tref = -0.0000000000008768554;
a6O2Tref = -1005.249;
if Tmp >= 1000
    a1O2T2 = 3.697578;
    a2O2T2 = 0.0006135197;
    a3O2T2 = -0.00000012588420;
    a4O2T2 = 0.00000000001775281;
    a5O2T2 = -0.0000000000000011364354;
    a6O2T2 = -1233.9301;
elseif Tmp < 1000
    a1O2T2 = 3.212936;
    a2O2T2 = 0.001127486;
    a3O2T2 = -0.0000005756150;
    a4O2T2 = 0.0000000013138773;
    a5O2T2 = -0.0000000000008768554;
    a6O2T2 = -1005.249;
end

% N2 %lec 4-6 slide 35
a1N2Tref = 3.298677;
a2N2Tref = 0.0014082404;
a3N2Tref = -0.000003963222;
a4N2Tref = 0.000000005641515;
a5N2Tref = -0.000000000002444854;
a6N2Tref = -1020.8999;

if Tmp >= 1000
    a1N2T2 = 2.926640;
    a2N2T2 = 0.0014879768;
    a3N2T2 = -0.0000005684760;
    a4N2T2 = 0.00000000010097038;
    a5N2T2 = -0.000000000000006753351;
    a6N2T2 = -922.7977;
elseif Tmp < 1000
    a1N2T2 = 3.298677;
    a2N2T2 = 0.0014082404;
    a3N2T2 = -0.000003963222;
    a4N2T2 = 0.000000005641515;
    a5N2T2 = -0.000000000002444854;
    a6N2T2 = -1020.8999;
end

% CO2
a1CO2Tref = 2.275724;
a2CO2Tref = 0.009922072;
a3CO2Tref = -0.000010409113;
a4CO2Tref = 0.000000006866686;
a5CO2Tref = -0.000000000002117280;
a6CO2Tref = -48373.14;
if Tmp < 1000
    a1CO2T2 = 2.275724;
    a2CO2T2 = 0.009922072;
    a3CO2T2 = -0.000010409113;
    a4CO2T2 = 0.000000006866686;
    a5CO2T2 = -0.00000000002117280;
    a6CO2T2 = -48373.14;
elseif Tmp >= 1000
    a1CO2T2 = 4.453623;
    a2CO2T2 = 0.003140168;
    a3CO2T2 = -0.0000012784105;
    a4CO2T2 = 0.0000000002393996;
    a5CO2T2 = -0.000000000000016690333;
    a6CO2T2 = -48966.96;
end

% H2O
a1H2OTref = 3.386842;
a2H2OTref = 0.003474982;
a3H2OTref = -0.000006354696;
a4H2OTref = 0.000000006968581;
a5H2OTref = -0.000000000002506588;
a6H2OTref = -30208.11;
if Tmp < 1000
    a1H2OT2 = 3.386842;
    a2H2OT2 = 0.003474982;
    a3H2OT2 = -0.000006354696;
    a4H2OT2 = 0.000000006968581;
    a5H2OT2 = -0.000000000002506588;
    a6H2OT2 = -30208.11;
elseif Tmp >= 1000
    a1H2OT2 = 2.672145;
    a2H2OT2 = 0.003056293;
    a3H2OT2 = -0.0000008730260;
    a4H2OT2 = 0.00000000012009964;
    a5H2OT2 = -0.000000000000006391618;
    a6H2OT2 = -29899.21;
end

% h0_fTref
h0_fuelTref = 4184*(a1*thetaTref + a2*(thetaTref^2)/2 + a3*(thetaTref^3)/3 + a4*(thetaTref^4)/4 - a5*thetaTref^-1 + a6);
h0_O2Tref = 8.314*(a1O2Tref*Tref + a2O2Tref*(Tref^2)/2 + a3O2Tref*(Tref^3)/3 + a4O2Tref*(Tref^4)/4 + a5O2Tref*(Tref^5)/5 + a6O2Tref);
h0_N2Tref = 8.314*(a1N2Tref*Tref + a2N2Tref*(Tref^2)/2 + a3N2Tref*(Tref^3)/3 + a4N2Tref*(Tref^4)/4 + a5N2Tref*(Tref^5)/5 + a6N2Tref);
h0_CO2Tref = 8.314*(a1CO2Tref*Tref + a2CO2Tref*(Tref^2)/2 + a3CO2Tref*(Tref^3)/3 + a4CO2Tref*(Tref^4)/4 + a5CO2Tref*(Tref^5)/5 + a6CO2Tref);
h0_H2OTref = 8.314*(a1H2OTref*Tref + a2H2OTref*(Tref^2)/2 + a3H2OTref*(Tref^3)/3 + a4H2OTref*(Tref^4)/4 + a5H2OTref*(Tref^5)/5 + a6H2OTref);

% h0_T1
h0_fuelT1 = 4184*(a1*thetaT1 + a2*(thetaT1^2)/2 + a3*(thetaT1^3)/3 + a4*(thetaT1^4)/4 - a5*thetaT1^-1 + a6);
h0_O2T1 = 8.314*(a1O2Tref*T + a2O2Tref*(T^2)/2 + a3O2Tref*(T^3)/3 + a4O2Tref*(T^4)/4 + a5O2Tref*(T^5)/5 + a6O2Tref);
h0_N2T1 = 8.314*(a1N2Tref*T + a2N2Tref*(T^2)/2 + a3N2Tref*(T^3)/3 + a4N2Tref*(T^4)/4 + a5N2Tref*(T^5)/5 + a6N2Tref);

% cp_T1
cp_fuelT1 = 4.184*(a1 + a2*thetaT1 + a3*thetaT1^2 + a4*thetaT1^3 + a5*thetaT1^-2);

% h0_T2 TmpER1
h0_CO2T4 = 8.314*(a1CO2T2*Tmp(i) + a2CO2T2*(Tmp(i)^2)/2 + a3CO2T2*(Tmp(i)^3)/3 + a4CO2T2*(Tmp(i)^4)/4 + a5CO2T2*(Tmp(i)^5)/5 + a6CO2T2);
h0_H2OT4 = 8.314*(a1H2OT2*Tmp(i) + a2H2OT2*(Tmp(i)^2)/2 + a3H2OT2*(Tmp(i)^3)/3 + a4H2OT2*(Tmp(i)^4)/4 + a5H2OT2*(Tmp(i)^5)/5 + a6H2OT2);
h0_O2T4 = 8.314*(a1O2T2*Tmp(i) + a2O2T2*(Tmp(i)^2)/2 + a3O2T2*(Tmp(i)^3)/3 + a4O2T2*(Tmp(i)^4)/4 + a5O2T2*(Tmp(i)^5)/5 + a6O2T2);
h0_N2T4 = 8.314*(a1N2T2*Tmp(i) + a2N2T2*(Tmp(i)^2)/2 + a3N2T2*(Tmp(i)^3)/3 + a4N2T2*(Tmp(i)^4)/4 + a5N2T2*(Tmp(i)^5)/5 + a6N2T2);

% h0_T2 TmpER8
%h0_CO2T2ER8 = 8.314*(a1CO2T2*TmpER8 + a2CO2T2*(TmpER8^2)/2 + a3CO2T2*(TmpER8^3)/3 + a4CO2T2*(TmpER8^4)/4 + a5CO2T2*(TmpER8^5)/5 + a6CO2T2);
%h0_H2OT2ER8 = 8.314*(a1H2OT2*TmpER8 + a2H2OT2*(TmpER8^2)/2 + a3H2OT2*(TmpER8^3)/3 + a4H2OT2*(TmpER8^4)/4 + a5H2OT2*(TmpER8^5)/5 + a6H2OT2);
%h0_O2T2ER8 = 8.314*(a1O2T2*TmpER8 + a2O2T2*(TmpER8^2)/2 + a3O2T2*(TmpER8^3)/3 + a4O2T2*(TmpER8^4)/4 + a5O2T2*(TmpER8^5)/5 + a6O2T2);
%h0_N2T2ER8 = 8.314*(a1N2T2*TmpER8 + a2N2T2*(TmpER8^2)/2 + a3N2T2*(TmpER8^3)/3 + a4N2T2*(TmpER8^4)/4 + a5N2T2*(TmpER8^5)/5 + a6N2T2);

end

