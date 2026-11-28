
clear all;          % clear workspace 
close all;          % close all figures 
clc;                % clear command window


%% Givens

epsilon = 16;       % compression ratio of cycle
eq = 0.8;           % equivalence ratio (fuel lean since less than 1)
T1 = 298;           % K, the temperature at point 1
p1 = 100;           % kPa, the pressure at point 1
n1 = 1.38;          % the polytropic index for compression (air)
n2 = 1.25;          % the polytropic index for expansion (combustion products)
R_u = 8.314;        % kj/(kmolK), universal gas constant 
 
In = 2;         % Initialize input 
while In ~=1 && In ~=0 % Begin Input While loop 
In = input('Would you like to keep the maximum pressure of 5500 kPa? type 1 if yes and 0 if no and hit enter: '); % Ask user input for keeping pressure value

if  In==1     % Begin if statement if user input was yes
pMAX = 5500;  % kPa, max pressure case 

elseif In==0  % Begin if statement if user input was no
In2 = input('What value would you like--please put in a real number value: ')
pMAX= In2;    % kPa, max pressure case
else          % Else statement
    fprintf('INPUT ERROR, PLEASE TRY AGAIN\n')% Display error

end           % end the if statement 
end           % end while loop 

fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n') % print line to command window 
fprintf('                                                                        Full Cycle Analysis\n')                                                                                                 % print text to command window
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n') % print line to command window
fprintf(' \n') % Print space to command window 
fprintf(' \n') % Print space to command window


% Define the fuel and oxidizer types
in3=0; % Initialize input 3
while in3 ~=1 && in3 ~=2 && in3 ~=3 % Begin Input while loop 
in3 = input('Which Fuel Should be Used? 1 for Gasoline 1(C8.26H15.5), 2 for Gasoline2(C7.76H13.1) or 3 for pure Hydrogen (H2): '); % Begin user input statement 
if in3 == 1             % Begin if statement if user input is 1
    fuel = 'Gasoline1'; % The fuel to be used
elseif in3 == 2         % Begin if statement if user input is 2
    fuel = 'Gasoline2'; % The fuel to be used
elseif in3 == 3         % Begin if statement if user input is 3
    fuel = 'Hydrogen';  % The fuel to be used
else                    % Else statement
fprintf('INPUT ERROR, PLEASE TRY AGAIN\n')% Display error
 
end                     % end if statement
end                     % end while loop




oxidizer = 'Air';       % Define the oxidizer used

if strcmp(fuel,'Gasoline1')     % Begin if statement if gasoline 1 is chosen
    MW_fuel = 114.8;            % kg/kmol, Define the molecular weight of the fuel 
elseif strcmp(fuel,'Gasoline2') % Begin if statement if gasoline 2 is chosen
    MW_fuel = 106.4;            % kg/kmol, Define the molecular weight of the fuel
elseif strcmp(fuel,'Hydrogen')  % Begin if statement if pure hydrogen is chosen
    MW_fuel = 2;                % kg/kmol, Define the molecular weight of the fuel
end                             % end the if statement 

%% Problem I.1


fprintf('\n ') % Print space to command window
fprintf('\n ') % Print space to command window


if eq<1 % Begin if statement when equivalence ratio is less than one
    fprintf('                                       This is a fuel lean process since the equivalence ratio is less than 1.\n')     % Print to the command window that the fuel is lean 
if eq>1 % Begin if statement when equivalence ratio is greater than one
    fprintf('                                       This is a fuel rich process since the equivalence ratio is greater than 1.\n')  % Print to the command window that the fuel is lean
if eq==1 % Begin if statement when equivalence ratio is equal to one
    fprintf('                                       This is a stoichiometric process since the equivalence ratio is equal to 1.\n') % Print to the command window that the fuel is lean
end % end first if statement
end % end second if statement
end % end third if statement



if strcmp(fuel,'Hydrogen') % Begin if statement if Hydrogen is chosen as the fuel type
    x=0;                   % Define the number of carbon atoms in the fuel
    y=2;                   % Define the number of hydrogen atoms in the fuel 
    astoich = (x+(1/4)*y); % Find the stoichiometric combustion constant
    a=astoich/eq;          % Find the coefficient for air and nitrogen 
    b =x;                  % Find the coefficient for CO2
    c =y/2;                % Find the coefficient for H2O
    d = a-b-(c/2);         % Find the coefficient for O2
    n = a*3.76;            % Find the value for the multiplication for the coefficient for air with 3.76
    a_mat = [a,b,c,d,n];   % Store all chemical equation values into matrix 
elseif strcmp(fuel,'Gasoline1') % Begin if statement for gasoline 1 fuel type 
    x=8.26;                     % Define the number of carbon atoms in the fuel 
    y=15.5;                     % Define the number of hydrogen in the fuel 
    [a,b,c,d,n] = Lean_Stoi_coeff(eq,fuel,'Air');  % Function that outputs the coefficents for 1 mol of fuel, matrix is [air,co2,h2o,o2,n2]
    a_mat = [a,b,c,d,n];        % Store all chemical equation values into matrix 
elseif strcmp(fuel,'Gasoline2') % Begin if statement for gasoline 2 fuel type 
    x=7.76;                     % Define the number of carbon atoms in the fuel 
    y=13.1;                     % Define the number of hydrogen atoms in the fuel 
    [a,b,c,d,n] = Lean_Stoi_coeff(eq,fuel,'Air');  % Function that outputs the coefficents for 1 mol of fuel, matrix is [air,co2,h2o,o2,n2]
    a_mat = [a,b,c,d,n];        % Store all chemical equation values into matrix
end                             % end if statement 
fprintf(' \n')                  % Create space in output 
fprintf('                                     C%4.2fH%4.2f + %4.2f(O2 + 3.76N2) -->  %4.2f(CO2) + %4.2f(H2O) + %4.2f(O2) + %4.2f(3.76N2)\n', x, y, a, b, c, d,a) % Output chemical equation to command window 


%% Problem I.2

% First determine the variables at certain ponts in the cycle
% Properties at point 1

R_air=8.313/28.96;         % kJ/kgK, ratio of specific heats for air (compression)
sv1=(R_air*T1)/(p1);       % (m^3)/kg, specific volume at point 1


% Properties at point 2

sv2=sv1/epsilon;           % (m^3)/kg, specific volume at point 2
p2=(epsilon^n1)*p1;        % kPa, pressure at point 2
T2=(epsilon^(n1-1))*T1 ;   % K, temperature at point 2 

[h_LHV_P,h_LHV_fuel] = h_LHV(eq,fuel,oxidizer,T2);  % h_LHV_P is the heating value per kg of products
                                                    % h_LHV_fuel is heating
                                                    % value per kg of fuel
LHV=h_LHV_fuel-h_LHV_P;                             % Final LHV value 


%% Problem I.3

T_test = EvaluateT(T2,h_LHV_P,'cv',a_mat);  % find the maximum temp if all of the heat goes into the constant volume heating

if (T_test/T2) <= pMAX/p2    % Begin if statement 
    alpha = T_test/T2;       % Find ratio of pressure increase value 
    T3 = T_test;             % K, Find temperature at point 3
    T4 = T_test;             % K, Find temperature at point 4 
    beta = 1;                % State beta value knowing previous assumption
elseif (T_test/T2) > pMAX/p2 % Begin if statement 
    alpha = pMAX/p2;         % Find the ratio of pressure increase value 
    n_mix = b + c + d + n;   % moles, Sum all the moles in chemical equation 
    
    MW_av = ((b*44 + c*18 + d*32 + n*28)/n_mix); % kg/kJ, Detemine the molecular weight to all the moles in chemical equation 
    
    T3 = T2*alpha; % K, Define Temperature at point 3 
  
    h_LHV_change = (( (b/n_mix) * ( (h_species('CO2',T3) - h_species('CO2',T2) )) + ...
    (c/n_mix)*( (h_species('H2O',T3) - h_species('H2O',T2)) ) + ...
    (d/n_mix)*( (h_species('O2',T3) - h_species('O2',T2)) ) + ...
    (n/n_mix)*( (h_species('N2',T3) - h_species('N2',T2)) )) - R_u*(T3 - T2))/MW_av; % The entalpy for the lower heating value 
    
    T4 = EvaluateT(T3,h_LHV_P - h_LHV_change ,'cp',a_mat);  % K, Temperature at point 4 
    beta = T4/T3; % Determine the value of the ratio of constant volume 
end % End the if statement 



%% Problem I.4


% Determine the rest of the variables at locations 4,5, and 6

% Properties at 3

sv3=sv1/epsilon ;                                    % (m^3)/kg, specific volume at point 3
p3 = alpha*((epsilon^n1)*p1)  ;                      % kPa, pressure at point 3  


% Properties at 4

sv4=(beta*sv1)/(epsilon);                            % (m^3)/(kg), specific volume at point 4 
p4 = alpha*((epsilon^n1)*p1);                        % kPa, pressure at point 4 

% Properties at 5

sv5=sv1;                                             % (m^3)/(kg), specific volume at point 5  
p5=((beta/epsilon)^(n2))*alpha*(epsilon^n1)*p1;      % kPa, pressure at point 5 
T5 = ((beta/epsilon)^(n1))*alpha*(epsilon^(n1))*T1;  % K, temperature at point 5  


% Find the Net work of the cycle using the work function 
Wcycle=Work(sv1,sv4, p1, epsilon, alpha, beta, n1, n2); % kJ/kg, define the work of the full cycle 




%% Problem I.5


% Find the mean effective pressure using the mep function
mep = Wcycle/(sv1-sv2);  % kPa, the mean effective pressure of the process 


%% Problem I.6


% Find the thermal efficiency of the cycle
Eta = Wcycle/h_LHV_P; % The thermal efficiency of the cycle

%% Problem I.7 


m_dot=55/(Eta*Wcycle);  % kg/s, the amount of fuel injected to attain a 55kW of power


%% Plot p-v diagram 
v_1_2 = (sv2:0.0005:sv1);               % m^3/kg, setup iteration for specific volume  
p_1_2 = zeros(1,length(v_1_2));         % kPa, setup iteration for pressure 
for i = 1:length(v_1_2)                 % Set up for statement  
    p_1_2(i) = p1* (sv1/v_1_2(i)) ^ n1; % kPa, pressure calculation 
end                                     % end the for loop 

v_4_5 = (sv4:0.0005:sv5);               % m^3/kg, setup iteration for specific volume 
p_4_5 = zeros(1,length(v_4_5));         % kPa, setup iteration for pressures 
for i = 1:length(v_4_5)                 % Begin for loop 
    p_4_5(i) = p4* (sv4/v_4_5(i)) ^ n2; % kPa, Pressure calculation 
end                                     % end the for loop
        

table([sv1,sv2,sv3,sv4,sv5],[p1,p2,p3,p4,p5]); % generate table of values
figure                                         % Initialize figure          
hold on                                        % Tell the figure to hold onto the previous inputs 
plot(v_1_2,p_1_2)                              % Plot the specific volume values to the pressures for points 1 to 2 
plot([sv2,sv3,sv4],[p2,p3,p4])                 % Plot the specific volume values to the pressures for points 2, to 3,to 4 
plot(v_4_5,p_4_5)                              % Plot the specific volume values to the pressures for points 4 to 5 
plot([sv5,sv1],[p5,p1])                        % Plot the specific volume values to the pressures for points 5 to 1
xlabel('Specific Volume, v [m^3/kg]');         % Label the x axis 
ylabel('Pressure, P [kPa]');                   % Label the y axis 
title('P vs v');                               % Title the figure 
grid on                                        % Generate grid on plot 


Prop = {'Lower Heating Value At T2';'Alpha';'Beta';'Net Work of the Cycle';'Mean Effective Pressure of the Cycle';'Thermal Efficiency of the Cycle';'Kg of Fuel per sec'}; % Generate table Variable names 
Value = [LHV,alpha,beta,Wcycle,mep,Eta,m_dot]';           % Store values in an array
Units = ["kJ/kg","NA","NA","kJ/Kg","kPa","NA","kg/sec"]'; % Store the unit names in an array 
T = table(Value,Units,'RowNames',Prop)                    % Generate a table of final values 
























%% FUNCTIONS

%% Lean stoichiometric coefficient

% The purpose of this function is to determine the chemical equation of the
% combustion process
function [a,b,c,d,n] = Lean_Stoi_coeff(eq_ratio,fuel,ox) 
%Takes a equiv ratio, fuel and oxidizer and balances chemical equation
%Valid fuels; 'Methane','Propane','Diesel','Gasoline1','Gasoline2'
%Valid Oxydizers; 'Air','Oxygen'
%
fuel_mat = [{'Methane'},1 ,4 ; {'Propane'}, 3, 8 ; {'Diesel'}, 10.8, 18.7; {'Gasoline1'}, 8.26, 15.5; {'Gasoline2'}, 7.76, 13.1]; % Define matrix of fuel types
for i = 1:length(fuel_mat(:,1))         % Begin for loop
    if strcmp(fuel,fuel_mat(i,1)) == 1  % Begin if statement fo gasoline 1
        x = cell2mat(fuel_mat(i,2));    % State value for number of carbon atoms
        y = cell2mat(fuel_mat(i,3));    % State value for number of hydrogen atoms
    end                                 % End if statement 
end                                     % End for loop 
if eq_ratio > 1                      % Begin if statement 
    error("Combustion is not Lean")  % If equivalence ratio is greater then 1, print error to command window  
end                                  % End if statement 
if strcmp(ox,'Air')            % Begin if statment if air is chosen 
    a_stoi = x +(y/4);         % Find stoichiomentric coefficient 
    a = a_stoi / eq_ratio;     % Find value for moles of air (also used in nitrogen calc in products)
    b = x;                     % Find value for moles of CO2
    c = y/2;                   % Find value for moles of H2O
    d = 0.5*(2*a - 2*b - c);   % Find value for moles of O2 
    n = a*3.76;                % Find number of moles for nitrogen  
elseif strcmp(ox,'Oxygen')   % Begin if statement if oxygen is used 
    a_stoi = x + (y/4);      % Determine stoichiometric coefficient  
    a = a_stoi / eq_ratio;   % Determine moles of air (also used in nitrogen calc in products)
    b = x;                   % Determine moles of CO2 
    c = y/2;                 % Determine moles of H2O
    d = 0.5*(2*a - 2*b - c); % Determine moles of O2 
    n = 0;                   % Find the number of moles for nitroge 
end                     % end if statement 
end                     % end for loop 



%% Enthalpy of fuel

% This function determines the enthalpy for a fuel (hydrocarbons only)
% Uses table for fuel given in lab document
% returns values in Kj/Kmol

function [h] = h_fuel(Sp,T)
%% For finding the enthalpy of Fuels only
% Uses table for fuel given in lab document
% returns values in Kj/Kmol

%% Find the correct a values for the required fuel
if strcmp('Methane',Sp) % Begin if statement if methane is chosen as fuel 
    a_mat = [-0.29149,26.327,-10.610,1.5656,0.16573,-18.331]; % Set up matrix of coefficients

elseif strcmp('Propane',Sp) % Begin if statement if Propane is chosen as fuel
    a_mat = [-1.4867,74.339,-39.065,8.0543,0.01219,-27.313]; % Set up matrix of coefficients

elseif strcmp('Diesel',Sp) % Begin if statement if diesel is chosen as fuel
    a_mat = [-9.1063,246.97,-143.74,32.329,0.0518,-50.128]; % Set up matrix of coefficients

elseif strcmp('Gasoline1',Sp) % Begin if statement if gasoline 1 is chosen as fuel
    a_mat = [-24.078,256.63,-201.68,64.750,0.5808,-27.562]; % Set up matrix of coefficients

elseif strcmp('Gasoline2',Sp) % Begin if statement if gasoline 2 is chosen as fuel
    a_mat = [-22.501,227.99,-177.26,56.048,0.4845,-17.578]; % Set up matrix of coefficients
end  % end if statement 


%% Calculate h with formula given on table
th = T/1000;    % Determine value for theta 
h = 4148*(a_mat(1)*th + a_mat(2)*th^2*(1/2) + a_mat(3)*th^3*(1/3) + a_mat(4)*th^4*(1/4) - a_mat(5)*th^-1 + a_mat(6)); % kJ of fuel/kg, Find the final enthalpy value 
end % end Function 




%% Evaluate T

% This function determines the alpha and beta values 

function Tf = EvaluateT(T0,h_LHV_P,cp_cv,coeff_mat) % Begin function 
 
b = coeff_mat(2);   % Determine coefficient 
c = coeff_mat(3);   % Determine coefficient 
d = coeff_mat(4);   % Determine coefficient 
n = coeff_mat(5);   % Determine coefficient


n_mix = b + c + d + n;  % find number of moles in full mixture 
MW_av = (b*44 + c*18 + d*32 + n*28)/n_mix; % Determine molecular weight of mixture 

if strcmp(cp_cv,'cv')           % Begin if statement if specific heat at constant volume 
    R_u = 8.314 ;               % kJ/kmolK, Output universal gas constant 
elseif strcmp(cp_cv,'cp')       % Begin if statement if specific heat at constant pressure 
    R_u = 0;                    % kJ/kmolK, Output zero
end                             % end if statement 
        
error = 1;          % Let the error be 1
T = T0 * 2;         % Let T be Tinitial times 2 
z = 1 ;
LHS = h_LHV_P * MW_av;  
RHS = (b/n_mix)*(h_species('CO2',T) - h_species('CO2',T0)) + ...
    (c/n_mix)*(h_species('H2O',T) - h_species('H2O',T0)) + ...
    (d/n_mix)*(h_species('O2',T) - h_species('O2',T0)) + ...
    (n/n_mix)*(h_species('N2',T) - h_species('N2',T0)) - R_u*(T - T0);

while abs(error) > 0.001 && z < 10000 % Begin While loop 
    
    error = RHS/LHS - 1 ; % State error bound
    tol_mat(z) = error;
       if error > 0       % Begin if statement  
           T = T - 10*abs(error);
           RHS = (b/n_mix)*((h_species('CO2',T) - h_species('CO2',T0))) + ...
    (c/n_mix)*((h_species('H2O',T) - h_species('H2O',T0))) + ...
    (d/n_mix)*((h_species('O2',T) - h_species('O2',T0))) + ...
    (n/n_mix)*((h_species('N2',T) - h_species('N2',T0))) - R_u*(T - T0);
           
       elseif error < 0  % Continue if statement  
           T = T + 10*abs(error);
           RHS = (b/n_mix)*((h_species('CO2',T) - h_species('CO2',T0))) + ...
    (c/n_mix)*((h_species('H2O',T) - h_species('H2O',T0))) + ...
    (d/n_mix)*((h_species('O2',T) - h_species('O2',T0))) + ...
    (n/n_mix)*((h_species('N2',T) - h_species('N2',T0))) - R_u*(T - T0);
           
       end              % end if statment 
       
    
    z = z + 1;          % Increase counting variable by 1
end                     % end while loop 

Tf = T;                 % State final temperature 
end                     % End Function 





%% Enthalpy of Species 

% This function determines the enthalpy for the species only
function [h] = h_species(Sp,T) % Begin function 
% For finding the enthalpy of given species only
% Uses table for species given in lab document
% returns values in Kj/Kmol

% Find the correct a values for the required fuel
if T > 5000 || T < 298 % Begin if statemetn 
     %error('temp out of valid range') % Print temp error 
end  % end if statement
if strcmp('CO2',Sp)     % Begin if statement for CO2
    if T > 1000         % Begin temperature range if statement 
        a_mat = [0.04453623E+02,0.03140168E-01,-0.12784105E-05,0.02393996E-08,-0.16690333E-13,-0.04896696E+06]; % Define constants 
    elseif T <= 1000    % Begin if statement for temperature range 
        a_mat = [0.02275724E+02,0.09922072E-01,-0.10409113E-04,0.06866686E-07,-0.02117280E-10,-0.04837314E+06]; % Define constants
    end                 % End if statement 
end                     % end if statement  
if strcmp('H2',Sp)      % Begin if statement if H2 is chosen 
    if T > 1000         % Begin if statemet temperature range 
        a_mat = [0.02991423E+02,0.07000644E-02,-0.05633828E-06,-0.09231578E-10,0.15827519E-14,-0.08350340E+04]; % Define constants
    elseif T <= 1000    % Begin if statement temperature range 
        a_mat = [0.03298124E+02,0.08249441E-02,-0.08143015E-05,-0.09475434E-09,0.04134872E-11,-0.10125209E+04]; % Define constants
    end                 % end if statement 
end                     % end if statement 
if strcmp('H2O',Sp)     % Begin if statement if H2O is chosen 
    if T > 1000         % Begin if statement for temperature range 
        a_mat = [0.02672145E+02, 0.03056293E-01 ,-0.08730260E-05 ,0.12009964E-09 ,-0.06391618E-13 ,-0.02989921E+06 ]; % Define constants 
    elseif T <= 1000    % Begin if statement 
        a_mat = [0.03386842E+02, 0.03474982E-01, -0.06354696E-04, 0.06968581E-07, -0.02506588E-10, -0.03020811E+06];  % Define constants 
    end                 % End if statement 
end                     % End if statement 
if strcmp('N2',Sp)      % Begin if statement for N2
    if T > 1000         % Begin if statement for temperature range 
        a_mat = [0.02926640E+02,0.14879768E-02,-0.05684760E-05,0.10097038E-09,-0.06753351E-13,-0.09227977E+04];  % Define constants
    elseif T <= 1000
        a_mat = [0.03298677E+02,0.14082404E-02,-0.03963222E-04,0.05641515E-07,-0.02444854E-10,-0.10208999E+04];  % Define constants 
    end                 % End if statement 
end                     % End if statement 
if strcmp('O2',Sp)      % Begin if statement for O2
    if T > 1000         % Begin if statement for temperature range 
        a_mat = [0.03697578E+02,0.06135197E-02,-0.12588420E-06,0.01775281E-09,-0.11364354E-14,-0.12339301E+04]; % Define constants 
    elseif T <= 1000    % Begin if statement 
        a_mat = [0.03212936E+02,0.11274864E-02,-0.05756150E-05,0.13138773E-08,-0.08768554E-11,-0.10052490E+04]; % Define constants
    end                 % End if statement 
end                     % End if statement 
%% Calculate h with formula given in class
Ru = 8.314; % kJ/kmolK, Define the universal gas constant 
h = Ru*(a_mat(6) + a_mat(1)*T + (1/2)*a_mat(2)*T^2 + (1/3)*a_mat(3)*T^3 + (1/4)*a_mat(4)*T^4 + (1/5)*a_mat(5)*T^5); % kJ/kg, determine the enthalpy non-fuel molecule 
end         % End function 






%% Lower Heating Value


function [h_LHV_P,h_LHV_fuel] = h_LHV(eq_ratio,fuel,ox,temp) % Begin function 
%   H_LHV_KGFUEL 
%   Fuel is STR, 'Gasoline1' or 'Gasoline2'
%   ox is oxidizer, 'Air' or 'Oxygen'
%   Takes state 2 information and calculates LHV based on 
%   1 KG of fuel and 1 KG of products
%
%
%   TODO:
%   Add MW calculator for xy fuel (or add to stoi calculator)
%   Check units

if strcmp(fuel,'Gasoline1')     % Begin if statement for gasoline 1
    MW_fuel = 114.8;            % kJ/kmol, Define molecular weight of fuel for gasoline 1
elseif strcmp(fuel,'Gasoline2') % Begin if statement for gasoline 2
    MW_fuel = 106.4;            % kJ/kmol, Define molecular weight of fuel for gasoline 2
elseif strcmp(fuel,'Hydrogen')  % Begin if statement for hydrogen as fuel 
    MW_fuel = 2;                % kJ/kmol, Define molecular weight of fuel for hydrogen 
end                             % End if statement 


if strcmp(fuel,'Hydrogen')      % Begin if statement if hydrogen 
    a = 0.625; b = 0 ; c = 1; d = 0.125; n = 2.35; % mole value matix for hydrogen fuel calulated by hand
else                            % else statement 
    [a,b,c,d,n]=  Lean_Stoi_coeff(eq_ratio,fuel,ox);  % Function that outputs the coefficents for 1 mol of fuel, matrix is [air,co2,h2o,o2,n2] 
end                             % End if statement

if strcmp(fuel,'Hydrogen')      % Begin if statement if hydrogen 
    HR = 1*h_species('H2',temp) + a*h_species('O2',temp); % kJ, Determine the enthalpy of reactants 
else        % Else statement 
    HR = 1*h_fuel(fuel,temp) + a*h_species('O2',temp);    % kJ, Determine the enthalpy of reactants 
end         % End if statement 
HP = b*h_species('CO2',temp) + c*h_species('H2O',temp) + d*h_species('O2',temp); % kJ, Determine the enthalpy of products 

% This is for 1 kmol of fuel

H = HR - HP; % kJ, Determine the difference between enthalpy of products and reactants

h_molar = H / 1;  % divide kJ by 1kmol of fuel

MW_P = b*44 + c*18 + d*32 + n*28 ;% kg/kmol Determine the molecular weight of the products

h_LHV_fuel = h_molar/MW_fuel ;    % kJ/kg of fuel, lower heating value (kJ per kg fuel)
h_LHV_P = h_LHV_fuel/( MW_P / MW_fuel); % kJ/kg, lower heating value (kJ per kg)

end     % End Function







%% Work

% The purpose of this function is to determine the total cycle work generated for a given fuel-air dual cycle

function [Wcycle] = Work(sv1,sv4, p1, epsilon, alpha, beta, n1, n2);   % Define function name and the required input and output variables 

delta = sv1/sv4;                                                       % Determine the expansion ratio of the fuel-air cycle

Wcycle= p1*sv1*(epsilon^(n1-1))*(alpha*(beta-1)+((alpha*beta)/(n2-1))*(1-(1/(delta^(n2-1))))-(1/(n1-1))*(1-(1/(epsilon^(n1-1))))); % kJ/kg, the total work generated for the air+fuel Cycle

end      % End function 

