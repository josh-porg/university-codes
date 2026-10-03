function h = enthalpycalculatorFuel(Species,Ti)

if strcmpi(Species,'Methane') == 1
   Fuel_aMatrix = [-0.29149,26.327,-10.61,1.5656,0.16573,-18.331];

elseif strcmpi(Species,'Propane') == 1
   Fuel_aMatrix = [-1.4867,74.339,-39.065,8.0543,0.01219,-27.313];

elseif strcmpi(Species,'Diesel') == 1
   Fuel_aMatrix = [-9.1063,246.97,-143.74,32.329,0.0518,-50.128];

elseif strcmpi(Species,'Hydrogen') ==1
   Fuel_aMatrix = [0.03298124*10^2,0.08249441*10^-2,-0.08143015*10^-5,-0.09475434*10^-9,0.04134872*10^-11,-0.10125209*10^4];
end

if strcmpi(Species,'Hydrogen') == 0
    theta = Ti/1000;
    Theta_Matrix = [theta,(1/2)*theta^2,(1/3)*theta^3,(1/4)*theta^4,-1*theta^-1,1];
    h = 4184*dot(Fuel_aMatrix,Theta_Matrix);
elseif strcmpi(Species,'Hydrogen') == 1
    Ti_Matrix = [Ti,(Ti^2)/2,(Ti^3)/3,(Ti^4)/4,(Ti^5)/5,1];
    h = 8.314*(dot(Ti_Matrix,Fuel_aMatrix));
end
