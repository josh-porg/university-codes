%AE571 Final Project Part II
%% 1

clear;
clc;

%Input Parameters
Fuel = 'Hydrogen';  %The fuel being used
ER = 0.8;           %The Equivilence Ratio
CompRatio = 16;     %The Engine Compression Ratio v1/v2
P1 = 100;           %Ambient Pressure in kPa
T1 = 298;           %Ambient Temperature in Kelvin
n1 = 1.38;          %Polytropic Index of Compression
n2 = 1.25;          %Polytropic Index of Expansion

%Function to calculate the chemical equation using the fuel and Equivilence
%ratio assuming fuel lean combustion
[a,b,c,d] = AtomBalanceLab1(Fuel,ER);
%CxHy + a*(O2 +3.76*N2) = b*CO2 + c*H2O + d*O2 + a*3.76*N2. Remember this
%formula every time you see me use a,b,c, or d. b is not used at all in
%the following calculations so if you do use a hydrocarbon and get a b, it
%won't break the code. It'll just give you a wrong answer.

%% 2

T2 = T1*CompRatio^(n1-1);   %Calculating T2 at end of compression in Kelvin
P2 = P1*CompRatio^n1;       %Pressure at point 2 in kPa

%Calculating enthalpy of all species at T2 in KJ/kmol.
hMolFuel = enthalpycalculatorFuel(Fuel,T2);
[hMolAir, hMolO2, hMolN2] = enthalpycalculatorAir(T2,T2);
    %This was made originally to give hMolAir at 1 temp and the others at a
    %different temp in Lab 1. Here everything is at T2 so both temps are T2
[hMolCO2, hMolH2O] = enthalpycalculatorProducts(T2);   %hMolCO2 isn't used in this 
%reaction but the calculator calculates it as the first output so we need
%to have it to get hH2O as an output. hMolCO2 will be ignored

hMolReactants = 1*hMolFuel + a*hMolAir;
hMolProducts = c*hMolH2O + d*hMolO2 + a*3.76*hMolN2;

%Lower Heating Value per kmol of Hydrogen in Kj/Kmol
hFuelMolLHV = hMolReactants - hMolProducts;

%Lower Heating Value per kg of Hydrogen by dividing hFuelMolLHV by the
%Molecular Weight of hydrogen gas
hFuelKgLHV = hFuelMolLHV/2;

%Calculating the mass of the products so I can convert hFuelKgLHV to be in
%terms of kg of Product
mProducts = c*18+d*32+a*3.76*28;

%Lower Heating Value of the reaction in Kj/Kg
hLHV = hFuelKgLHV*2/mProducts;

%% 3

                        %Extra Credit

EngineMax = 5.5;     %Max Sustainable Engine Pressure in MPa.

%Need to find the maximum temperature that this combustion can produce.
%i.e. the adiabatic flame temperature. I'm recycling some code for lab 1
%for this.

    error = (hMolProducts - hMolReactants)/hMolReactants*100;
    T = 1000;   %Initial Guess Temperature in Kelvin
    i = 1;
    while abs(error) > 0.1 %0.1 percent tolerance
    %as ER gets smaller, the iteration becomes extremely sensetive to
    %changes in temp, which error makes large changes to. This creates a
    %feedback loop that results in temps in the 10^150 range.
    %1000*ER makes the maximum acceptable error to add or subtract
    %smaller as it becomes more sensitive to it.
        if abs(error) < 1000
            if error > 0
                T = T + error;  
            elseif error < 0
                T = T - error;
            end
        elseif abs(error) > 1000 
            %If error becomes too large, it typically 
      % means the error is varying widly and I should go by a smaller
      % change in T. Exp term gets smaller as i increases
            if error > 0
                T = T + 0.5*T*ER*exp(-0.1*i);
            elseif error < 0
                T = T - 0.5*T*ER*exp(-0.1*i);
            end
        end

    [h_air,h_O2,h_N2] = enthalpycalculatorAir(298,T);
    [h_CO2,h_H2O] = enthalpycalculatorProducts(T);
    h_products = b*h_CO2 + c*h_H2O + d*h_O2 + a*3.76*h_N2;
    error = (h_products - hMolReactants)/hMolReactants*100;
%Counting Iterations    
    i = i + 1;
%Logging values per Iteration.
    IterationGraph(i) = i;
    TemperatureGraph(i) = T;
    ErrorGraph(i) = error;
%These two are just for debugging purposes
    ProductGraph(i) = h_products;
    ReactantGraph(i) = hMolReactants;
    
 %Crash and wrong answer Prevention
    if i >= 500
        break
    elseif T < 0
        break
    end
    end
%These are also for debugging
%subplot(1,2,1)
%plot(IterationGraph,TemperatureGraph)
%subplot(1,2,2)
%plot(IterationGraph,ErrorGraph)
MaxT = T;

%End of recycled Lab 1 code

%The maximum possible temperature is given by MaxT. Corresponding to this
%is a maximum possible pressure MaxP.  If the Engine Max Pressure EngineMax is 
%less than the MaxP, then the excess energy will be transfered into
%specific volume increase beta(v4/v3) since alpha and beta are inversely
%porportional. However, this means if EngineMax is larger than MaxP, all 
%the energy becomes pressure increase with no specific volume
%increase: i.e It becomes an Otto cycle with beta = 1. 
%So the code needs to check if EngineMax is greater than P3 if beta = 1.
%MaxP = R*T3/v3 = Rproducts*MaxT/v2 if beta = 1.

Ru = 8.3145; %Universal Gas Constant in Kj/KmolK

% Find Rreact (Gas constant of Reactants)
MWair = (32+3.76*28)/(1+3.76);  %air is already a mixture so it has to be
Mreact = (2+a*MWair);           %treated as a single species when finding 
MWreact = Mreact/(1+a);         %the Molecular Weight of the reactants
Rreact = Ru/MWreact;

% Use ideal gas law to find V1
V1 = (Rreact*T1)/P1;

% Find V2
V2 = V1/CompRatio;

%MW,mix = Sum(n,i*MW,i)/Sum(n,i).
MWmix = mProducts/(c+d+a*3.76);
Rmix = Ru/MWmix;

MaxP = Rmix*MaxT/V2;  %Maximum Possible Pressure due to the chemical
                      %reaction in kPa

if MaxP > EngineMax*10^3
    P3 = EngineMax*10^3;  %Converting EngineMax to kPa
    alpha = P3/P2;  %This is P3 at the max engine pressure
    T3 = alpha*T2;  %Temperature at Point 3 in Kelvin

    %beta = v4/v3.  v4 = R*T4/P4. R is Rmix, T4 = MaxT. P4 = P3. v3 = v2.
    %finding v4
    V4 = Rmix*MaxT/P3;
    
    %Finding Beta
    beta = V4/V2;

elseif MaxP < EngineMax*10^3 %if the maximum engine pressure is greater then it 
                        %becomes an Otto cycle
    P3 = MaxP;
    beta = 1;
    alpha = P3/P2;   %This is P3 at the max chemically possible pressure
    T3 = alpha*T2;
end






%% INPUT CHECK SECTION (2.4-2.7)
% use input alpha
% use input beta
% use input a from chemical equation
% given P1, T1
% given n1, n2
% given CompRatio epsilon

% Find V4
V4 = beta*V1/CompRatio;

% Solve for expansion ratio delta delta = (V5/V4) = (V1/V4)
delta = V1/V4;
sigma = 1;                                                          % set sigma to 1 for an ideal cycle assumption




%% 4. Calculate Net Work of Cycle
% Wcycle = sigma*P1*V1*(CompRatio^(n1-1))*(  alpha*(beta-1) + (alpha*beta/(n2-1))*(1 - 1/(delta^(n2-1)))- (1/(n1-1))*(1-(1/(CompRatio^(n1-1))))  )
% Use V1 and delta from previous sections

Coef_Wcycle = sigma*P1*V1*(CompRatio^(n1-1));                       % Calculatates the product of all the values in the coefficient for cycle Net Work CoefWcycle = sigma*P1*V1*(CompRatio^(n1-1))
Term1 = (alpha*(beta-1));                                           % 1st term in above equation
Term2 = (alpha*beta)/(n2-1);                                        % 2nd term in above equation
Term3 = (1-(1/(delta^(n2-1))));                                     % 3rd term in above equation
Term4 = (1/(n1-1));                                                 % 4th term in above equation
Term5 = (1 - (1/(CompRatio^(n1-1))));                                   % 5th term in above equation

Wcycle = Coef_Wcycle*(Term1 + (Term2*Term3) - (Term4*Term5));       % Total Net work of cycle (kJ)




%% 5. Calculate MEP of Cycle
% MEP = Wcycle/(Vmax-Vmin)    where Vmax = V1 and Vmin = V2
% Use Wcycle and V1 & V2 from previous sections

MEP = Wcycle/(V1-V2);                                               % Mean Effective Pressure (kPa) of the cycle


%% 6. Calculate Thermal Efficiency of Cycle
% Thermal Efficiency: eta_th = Wcycle/qin = Wcycle/hLHV
% Use Wcycle and hFuelKgLHV from previous sections


ThermalEfficiency = Wcycle/hLHV;




%% 7. Calculates necessary kg per second of fuel to produce 55 kJ per second
% CyclePower = WdotCycle = Wcycle/second
% Using given TARGET WdotCycle of 55 kW and previously calculated Thermal
% Efficiency, use the found TARGET heat in (qdot) and the previously
% calculated hFuelKgLHV and find the REQUIRED fuel input rate (mdotFuelin)

qdot = 55/ThermalEfficiency;                                        % qdot = WdotCycle/ThermalEfficiency   where WdotCycle = 55 kJ   (kJ/s)
mdotFuelin = qdot/hFuelKgLHV;                                       % mdotFuelin = qdot / hFuelKgLHV   (kg/s)



%% OUTPUT SECTION

ChemEQ = ['Chemical Equation:   H2 + ',num2str(a),'(O2 + 3.76*N2) <===> ',num2str(c),'*H2O + ',num2str(d),'*O2 + ',num2str(a*3.76),'*N2'];
disp(ChemEQ)

Name = {'LHV';'Alpha';'Beta';'Net Work';'MEP';'Thermal Efficiency';'Req. Fuel Flow Rate'};
Units = {'kJ/kg';'---';'---';'kJ/kg';'kPa';'---';'kg/s'};
Values = [hLHV;alpha;beta;Wcycle;MEP;ThermalEfficiency;mdotFuelin];
Output = table(Name,Values,Units)