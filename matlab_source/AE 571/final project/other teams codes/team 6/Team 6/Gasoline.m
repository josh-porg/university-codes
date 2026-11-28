%% GROUP 6 FINAL PROJECT - GASOLINE CODE
clear, clc
Ru = 8.3145;
Tref = 298;
%% DESIGN PARAMETERS
    % Instructions - Input Design Parameters. For testing, input parameters
    % directly

disp(['Considering a Four-Stroke Compression-Ignition (CI) Dual-Cycle' ...
    ' Reciprocating Engine with the Following Design Parameters:'])
disp(' ')
disp('CYCLE DESIGN PARAMETERS')

% epsilon = input('Compression Ratio: ');
       epsilon = 16;
       fprintf('Epsilon: %3.2f \n', epsilon)
% ER = input('Equivalence Ratio: ');
       ER = 0.8;
       fprintf('Equivalence Ratio: %3.2f \n', ER)
% T1 = input('Ambient Temperature (T1), K: ');
       T1 = 298;
       fprintf('Temperature @ State 1: %3.2f K \n', T1)
% P1 = input('Ambient Pressure (P1), Pa: ');
       P1 = 100E3;
       fprintf('Pressure @ State 1: %3.2f Pa \n', P1)
% n1 = input('Polytropic Index for Compression: ');
       n1 = 1.38;
       fprintf('Polytropic Index for Compression: %3.2f \n', n1)
% n2 = input('Polytropic Index for Expansion: ');
       n2 = 1.25;
       fprintf('Polytropic Index for Expansion: %3.2f \n', n2)

%% CHEMICAL REACTION EQUATION

    % FUEL ----------------------------------------------------------------
    % CxHy - Gasoline
      x = 8.26;
      y = 15.5;

    disp(' ')
    fprintf('Using Gasoline (C%3.2fH%3.1f): \n', x, y)

    % ATOM BALANCE FUNCTION -----------------------------------------------
    [a, b, c, d, e] = atombalance(x,y,ER);

    % CHEMICAL REACTION EQUATION ------------------------------------------
    if ER <= 1
        disp(' ')
        disp('CHEMICAL EQUATION: FUEL LEAN')
        fprintf('C%1.2fH%1.2f + %1.2f*(O2 + 3.76*N2) -> %1.2f*CO2 + %1.2f*H2O + %1.2f*O2 + %1.2f*(3.76*N2) \n', x,y,a,b,c,d,a)
    
    elseif ER >=1
        disp(' ')
        disp('CHEMICAL EQUATION: FUEL RICH')
        fprintf('C%1.2fH%1.2f + %1.2f*(O2 + 3.76*N2) -> %1.2f*CO2 + %1.2f*H2O + %1.2f*O2 + %1.2f*CO + %1.2f*(3.76*N2) \n', x,y,a,b,c,d,e,a)
    end

%% LOWER HEATING VALUE @T2

    % T2 VALUE ------------------------------------------------------------
    T2 = T1*epsilon^(n1-1);

    % REACTANTS -----------------------------------------------------------
    % Fuel
    [MW,h,cp] = propertycalculator('C8.26H15.5',T2);
    FuelMW = MW;                            % Molecular Weight of Substance, kg/kmol
    FuelCP = double(1*int(cp, Tref, T2));   % Molar cp Value of Substance, kJ
    Fuelcp_function = 1*cp;                 % Molar cp Function, kJ/K
    Fuel_h_Tref = 1*subs(h, Tref);          % Molar Enthalpy Value @298K of Substance, kJ
    Fuel_h = 1*h;                           % Molar Enthalpy Function of Substance, kJ

    % O2
    [MW,h,cp] = propertycalculator('O2',T2);
    O2RMW = MW;
    O2RCP = double(a*int(cp, Tref, T2));
    O2Rcp_function = a*cp;
    O2R_h_Tref = a*subs(h, Tref);
    O2R_h = a*h;

    % N2
    [MW,h,cp] = propertycalculator('N2',T2);
    N2RMW = MW;
    N2RCP = double(a*3.76*int(cp, Tref, T2));
    N2Rcp_function = a*3.76*cp;
    N2R_h_Tref = a*3.76*subs(h, Tref);
    N2R_h = a*3.76*h;

    % PRODUCTS ----------------------------------------------------------------
    % CO2
    [MW,h,cp] = propertycalculator('CO2',T2);
    CO2MW = MW;
    CO2CP = double(b*int(cp, Tref, T2));
    CO2cp_function = cp; %b
    CO2_h_Tref = b*subs(h, Tref);
    CO2_h = b*h;

    % H2O
    [MW,h,cp] = propertycalculator('H2O',T2);
    H2OMW = MW;
    H2OCP = double(c*int(cp, Tref, T2));
    H2Ocp_function = cp; %c
    H2O_h_Tref = c*subs(h, Tref);
    H2O_h = c*h;

    % O2
    [MW,h,cp] = propertycalculator('O2',T2);
    O2PMW = MW;
    O2PCP = double(d*int(cp, Tref, T2));
    O2Pcp_function = cp; %d
    O2P_h_Tref = d*subs(h, Tref);
    O2P_h = d*h;

    % CO
    [MW,h,cp] = propertycalculator('CO',T2);
    COMW = MW;
    COCP = double(e*int(cp, Tref, T2));
    COcp_function = cp; %e
    CO_h_Tref = e*subs(h, Tref);
    CO_h = e*h;

    % N2
    [MW,h,cp] = propertycalculator('N2',T2);
    N2PMW = MW;
    N2PCP = double(a*3.76*int(cp, Tref, T2));
    N2Pcp_function = cp; %3.76*a
    N2P_h_Tref = a*3.76*subs(h, Tref);
    N2P_h = a*3.76*h;

% MW ----------------------------------------------------------------------
sum_MW = (CO2MW + H2OMW + O2PMW + COMW + N2PMW); % Summation of MW, kg/kmol
sum_kg = (b*CO2MW + c*H2OMW + d*O2PMW + e*COMW + a*3.76*N2PMW); % Summation of Mass, kg
sum_Moles = (b + c + d + e + a*3.76); % Summation of Moles, kmol

% HR & HP -----------------------------------------------------------------
HR = (Fuel_h_Tref + O2R_h_Tref + N2R_h_Tref) + (FuelCP + O2RCP + N2RCP); % kJ
HP = (CO2_h_Tref + H2O_h_Tref + O2P_h_Tref + CO_h_Tref + N2P_h_Tref) + (CO2CP + O2PCP + COCP + N2PCP); % kJ
DeltaH = HR - HP; % kJ

% MOLAR HEATING VALUE -----------------------------------------------------
hcFuel_bar = DeltaH/1; % nFuel = 1; kJ/(kmol of fuel)

% HEATING VALUE OF FUEL ---------------------------------------------------
hcFuel = hcFuel_bar/(FuelMW); %kJ/(kg of fuel)
disp(' ')
disp('LOWER HEATING VALUES')
fprintf('Lower Heating Value of the Fuel: %3.2f kJ/kg_fuel \n', hcFuel)

% HEATING VALUE OF MIXTURE ------------------------------------------------
hcMix = hcFuel*(FuelMW/(sum_kg)); % kJ/kg
fprintf('Lower Heating Value of the Mixture: %3.2f kJ/kg_mix \n', hcMix)

% TOTAL HEAT FOR COMBUSTION -----------------------------------------------
hc = hcMix*sum_kg;

%% RATIO OF PRESSURE INCREASE AND VOLUME INCREASE

    % AlPHA ---------------------------------------------------------------
        % P2 CALCULATION
        P2 = P1*epsilon^(n1);
        % P3 CALCULATION
        beep
        disp(' ')
        P3 = input('Input MAX Engine Pressure, Pa -- (Default: 5.5E6 Pa): ');
        % P3 = 5.5E6;

        % P3 CALCULATION CHECK
        if P3 < P2
            beep
            disp(' ')
            disp('- ERROR: P3 Must be >= P2')
            fprintf('- P2 = %3.2f Pa \n', P2)
            return
        end
        

    % ALPHA CALCULATION
    alpha = P3/P2;

    % REACTANTS WITH ALPHA CORRECTION -------------------------------------
    % Since alpha*T2 may be higher than 1000 K, the coefficients 
    % in the enthalpy function will be different. To compensate 
    % for this difference, the following code will reanalyze the enthalpy 
    % function with the different coefficients. 

    % FUEL
    [~,h,~] = propertycalculator('C8.26H15.5',1000);
    Fuel_h_correction = 1*h;

    % O2R
    [~,h,~] = propertycalculator('O2',1000);
    O2R_h_correction = a*h;

    % N2R
    [~,h,~] = propertycalculator('N2',1000);
    N2R_h_correction = a*3.76*h;

    % PRODUCTS WITH ALPHA CORRECTION --------------------------------------
    % CO2
    [~,h,~] = propertycalculator('CO2',1000);
    CO2_h_correction = b*h; 

    % H2O
    [~,h,~] = propertycalculator('H2O',1000);
    H2O_h_correction = c*h; 

    % O2
    [~,h,~] = propertycalculator('O2',1000);
    O2P_h_correction = d*h; 

    % CO
    [~,h,~] = propertycalculator('CO',1000);
    CO_h_correction = e*h; 

    % N2
    [~,h,~] = propertycalculator('N2',1000);
    N2P_h_correction = a*3.76*h; 

    % BETA ----------------------------------------------------------------
        %  ENTHALPY OF REACTANTS & PRODUCTS T < 1000
        h_Products = (CO2_h + H2O_h + O2P_h + CO_h  + N2P_h); % kJ
        h_Reactants = (Fuel_h + O2R_h + N2R_h); % kJ

        % ENTHALPY OF REACTANTS & PRODUCTS WITH ALPHA CORRECTION T > 1000
        h_Products_correction = (CO2_h_correction + H2O_h_correction + O2P_h_correction + CO_h_correction + N2P_h_correction);
        h_Reactants_correction = (Fuel_h_correction + O2R_h_correction + N2R_h_correction);
     
        % RMIX OF REACTANTS
        yFuel = 1/sum_Moles;
        yO2R = a/sum_Moles;
        yN2R = (a*3.76)/sum_Moles;
        Molar_Mix = yFuel*FuelMW + yO2R*O2RMW + yN2R*N2RMW; % kg/kmol
        Rmix = Ru/Molar_Mix; % kJ/(kg*K)

        % qIN
        syms beta
        if alpha*T2 < 1000
            q1 = (subs(h_Reactants, alpha*T2) - subs(h_Reactants, T2));
        elseif alpha*T2 > 1000
            q1 = (subs(h_Reactants_correction, alpha*T2) - subs(h_Reactants, T2));
        end

        if alpha*T2 < 1000
            q2 = (subs(h_Products, alpha*beta*T2) - subs(h_Products, alpha*T2));
        elseif alpha*T2 > 1000
            q2 = (subs(h_Products_correction, alpha*beta*T2) - subs(h_Products_correction, alpha*T2));
        end

        % qIN BETA CALCULATION
        qin(beta) = q1 + q2 - Rmix*(alpha*T2 - T2)*sum_kg;

    % INITIAL VARIABLES ---------------------------------------------------
    qin_subbed = zeros(426,1);
    err = zeros(426,1);
    err(1,1) = 1;
    k = 1E-6;
    t = 1;
    tol = 0.001;
    disp('Calculating...')

    % ERROR LOOP FOR BETA -------------------------------------------------
    while err(t,1) > tol
        % Error Check
        qin_subbed(t,1) = subs(qin, k); % Substitutes Beta for k
        err(t+1,1) = hc - qin_subbed(t,1); % Checks Difference Between hc and qin

            if err(t+1,1) > 10000
                k = k + 1E-2;
            end
            if err(t,1) < 10000 && err(t,1) > 1000
                k = k + 1E-3;
            end
            if err(t+1,1) < 1000 && err(t,1) > 100
                k = k + 1E-5;
            end
            if err(t+1,1) < 100 && err(t,1) > 0.01
                k = k + 1E-9;
            end
            if err(t+1,1) < 0.01 && err(t,1) > 0
                k = k + 1E-10;
            end

            t = t+1;
                if t == 1E3
                    break
                end
    end

    % BETA CALCULATION ----------------------------------------------------
    beta = round(k,2);

    % CHECK ---------------------------------------------------------------
        % Instructions - Check Intercept of Both Plots
        figure('Name',"Comparison Between the Cycle's Heat of Combustion and Qin",'NumberTitle','off','Renderer', 'painters', 'Position', [100 100 1500 700]);
        fplot(qin, 'b')
        hold on
        fplot(hc, 'r'), title('hC & qIN Over Beta'), xlabel('beta (\beta)'), ylabel('Energy, (kJ)'), legend('qin', 'hc'), xlim([0,5]), grid on

% ALPHA AND BETA ----------------------------------------------------------
disp(' ')
disp('PV DIAGRAM VARIABLES')
fprintf('Alpha: %3.2f \n', alpha)
fprintf('Beta: %3.2f \n', beta)
%% NET WORK OF THE CYCLE

% Must First Determine V1 & Delta
    % V1 CALCULATION ------------------------------------------------------
    Rair = Ru/(28.9647*1E-3);
    V1 = (Rair*T1)/P1;

    % DELTA CALCULATION ---------------------------------------------------
    delta = epsilon/beta;
    
% WORK CALCULATION --------------------------------------------------------
W = (P1*V1*epsilon^(n1 - 1)*( alpha*(beta - 1) + ((alpha*beta)/(n2 - 1))*(1 - (delta^(n2 - 1))^(-1)) - ((n1 - 1)^(-1))*(1 - (epsilon^(n1 - 1))^(-1)) ))*1E-3; %kJ/kg

% WORK RELATIONSHIP
syms beta
Work(beta) = (P1*V1*epsilon^(n1 - 1)*( alpha*(beta - 1) + ((alpha*beta)/(n2 - 1))*(1 - ((epsilon/beta)^(n2 - 1))^(-1)) - ((n1 - 1)^(-1))*(1 - (epsilon^(n1 - 1))^(-1)) ))*1E-3; %kJ/kg
        figure('Name',"Comparison Between the Cycle's Qin & Work",'NumberTitle','off','Renderer', 'painters', 'Position', [100 100 1500 700]);
        fplot(qin, [0, 10], 'b')
        hold on
        fplot(Work, [0, 10], 'r'), title('Work & qIN Over Beta'), xlabel('Beta, (\beta)'), ylabel('Energy, (kJ)'), legend('qin','Work'), grid on

% NET WORK OF THE CYCLE ---------------------------------------------------
disp(' ')
disp('NET WORK OF THE CYCLE')
fprintf('Work: %3.2f kJ/kg \n', W)
%% MEP OF THE CYCLE

% Must First Determine V2, (Vmax)
    % V2 CALCULATION ------------------------------------------------------
    V2 = ((P1/P2)^(1/n1))*V1; % m^3/kg

% MEP CALCULATION ---------------------------------------------------------
mep = W/(V1-V2); % kPa

% MEP OF THE CYCLE --------------------------------------------------------
disp(' ')
disp('MEP OF THE CYCLE')
fprintf('mep: %3.2f kPa \n', mep)
%% THERMAL EFFICIENCY OF THE CYCLE

% EFFECIENCY CALCULATION --------------------------------------------------
n = (W/hcMix)*100;

% THERMAL EFFICIENCY OF THE CYCLE -----------------------------------------
disp(' ')
disp('THERMAL EFFICIENCY OF THE CYCLE')
fprintf('n: %3.2f Percent \n', n)
%% POWER

% Net Work of the Cycle = kJ/kg of Fuel
% Current Power = (Net Work of the Cycle) / (1 Sec)
% Desired Power = 55 kW

% NEEDED MASS FOR DESIRED WORK --------------------------------------------
neededmass = 55/W;

disp(' ')
disp('REQUIRED KG OF FUEL PER SEC TO PRODUCE 55KW')
fprintf('Required Mass Flow Rate: %3.3f kg/s \n', neededmass)

disp(' ')
disp('COMPLETE')
load gong
sound (y, Fs)

disp(' '), disp(' ')
disp('Code Created By: Oliver Gonzalez & Jason Oswald')