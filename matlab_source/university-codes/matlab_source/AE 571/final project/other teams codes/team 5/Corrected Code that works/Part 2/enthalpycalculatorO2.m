function h_O2i = enthalpycalculatorO2(Ti)
%Finding enthalpy of pure oxygen at a temperature from ~300K to 5000K
if Ti > 1000
   O2_Matrix = [0.03697578*10^2,0.06135197*10^-2,-0.1258842*10^-6,0.01775281*10^-9,-0.11364354*10^-14,-0.12339301*10^4];

elseif Ti <= 1000
    O2_Matrix = [0.03212936*10^2,0.11274864*10^-2,-0.0575615*10^-5,0.13138773*10^-8,-0.08768554*10^-11,-0.1005249*10^4];
end

Ti_Matrix = [Ti,(Ti^2)/2,(Ti^3)/3,(Ti^4)/4,(Ti^5)/5,1];
h_O2i = 8.314*(dot(Ti_Matrix,O2_Matrix));
end