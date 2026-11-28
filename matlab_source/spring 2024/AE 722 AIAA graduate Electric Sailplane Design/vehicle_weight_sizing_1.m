clc
clear all

% Aircraft Weight Sizing Using Dr. Roskam's Method

% Mission Profile 1: Launch Hermes at 65,000 ft
W_frac_hermes = 0.25;

W_toguess1 = 140000; % Guess take-off weight of aircraft


%W_pltot = (30000+30000*W_frac_hermes)*2; % Total payload is 2 Hermes rockets, each containing 30,000 of fuel and the structure weighs 20% of the fuel
W_pltot= (30000*W_frac_hermes)
W_crew = 175*4; % 4 crew required, each weighing 175 lbf

W_pl_crew = W_pltot+W_crew;

% Weight Fractions and other calculations to find fuel weight (W_f)
W_1_to = 0.99; % Estimated engine start W frac. from Roskam
W_2_1 = 0.99; % Estimated taxi W frac. from Roskam
W_3_2 = 0.995; % Estimated take-off W frac. from Roskam
W_4_3 = 0.98; % Estimated climb to 65,000 ft W frac. from Roskam

% Find weight after payload launch to use in Breguet Endurance Equation
W_4 = W_1_to*W_2_1*W_3_2*W_4_3*W_toguess1; % Climb weight (lbf) with payload
W_45 = W_4-W_pltot; % Weight (lbf) immediately after climb and payload launch

% Fuel weight used until payload launch
M_ff1 = W_1_to*W_2_1*W_3_2*W_4_3; % Mission fuel fraction
W_f_used1 = (1-M_ff1)*W_toguess1;

% Breguet Endurance Equation for loitering 8 min after payload launch
E_1min = 8; % Loiter for 8 minutes for Hermes to complete mission and be guided by crew
E_1 = E_1min/60; % Convert minutes to hours
L_D = 26; % L/D 26:1 assuming nlf P180 configuration
c_j = 0.34; % TSFC (lbf/lbf hr) assuming better than current engines (0.35) but using greener fuels so efficiency is not super improved

W_5 = W_45/(exp((E_1*c_j)/L_D)); % Loiter weight (lbf)
W_5_45 = W_5/W_45; % Estimated loitering W frac.

% W frac. continued
W_6_5 = 0.995; % Estimated descent W frac. +0.005 for efficiency from Roskam
W_7_6 = 1; % Estimated botched landing W frac. from Dr. Barrett
W_8_7 = 0.995; % Estimated small climb W frac. from Dr. Barrett

% Breguet Range Equation for divert to alternate 100 nmi away
R = 100; % 100 nmi range
M = 0.75; % Cruise speed as specified by objective function
a_fts = 1036.8; % 20,000 speed of sound ft/s (short climb altitude)
C_fts_kts = 3600/6076; % 3600 s/hr / 6076 ft/nmi
a = a_fts*C_fts_kts; % Convert speed of sound from ft/s to kts
V = M*a; % Divert to alternate velocity (kts)

W_9_8 = 1/exp((R*c_j)/(V*L_D)); % Rearranged B. Range Eqn. to find divert to alternate W frac.

% Breguet Endurance Equation for 45 min loiter after diversion per FAR 25
E_2min = 45; % 45 min loiter endurance
E_2 = E_2min/60; % Convert to hr
W_10_9 = 1/exp((E_2*c_j)/(L_D)); % Estimated loiter W frac.

% W frac. continued
W_11_10 = 0.995; % Estimated descent W frac. +0.005 for efficiency from Roskam and Dr. Barrett
W_12_11 = 1; % Estimated landing W frac. from Dr. Barrett
W_13_12 = 0.995; % Estimated taxi and shutdown W frac. +0.003 for efficiency form Roskam and Dr. Barrett

% Mission Fuel Fraction Equation to find weight of fuel used after payload launch
M_ff2 = W_5_45*W_6_5*W_7_6*W_8_7*W_9_8*W_10_9*W_11_10*W_12_11*W_13_12;
W_f_used2 = (1-M_ff2)*W_45; % Weight (lbf) of fuel used after payload launch
W_f_used = W_f_used1 + W_f_used2; % Total fuel used

% Tentative Operating Empty Weight (lbf) Equation
W_oe_tent = W_toguess1 - W_f_used - W_pl_crew

% Tentative Empty Weight (lbf) Equation and total fuel weight calculation
W_tfo = 0.01*W_f_used % Estimated trapped fuel and oil weight (lbf)
W_e_tent = W_oe_tent - W_tfo % Tentative empty weight (lbf)
W_f = W_f_used + W_tfo

W_4_5 = W_45 / W_4; % weight loss due to payload deploy
M_ff_landing = W_1_to*W_2_1*W_3_2*W_4_3* W_4_5 *W_5_45*W_6_5*W_7_6*W_8_7*W_9_8*W_10_9*W_11_10*W_12_11 % this is wieght fraction of everying lost during flight including fuel and hermes.

%% Mission 2 Profile: 3000 nmi Ferry
clc

W_frac_hermes = 0.25;
W_propellant = 30000*2;
W_pltot = W_propellant+(30000*W_frac_hermes)*2; % Total payload is 2 Hermes rockets, each weighing 333333 lbf
W_crew = 175*4; % 4 crew required, each weighing 175 lbf
W_toguess1 = 130000-W_pltot; %130000-W_pltot  Guess take-off weight of aircraft, assuming payload shells are being carried

% Weight Fractions and other calculations to find fuel weight (W_f)
W_1_to = 0.99; % Estimated engine start W frac. from Roskam
W_2_1 = 0.99; % Estimated taxi W frac. from Roskam
W_3_2 = 0.995; % Estimated take-off W frac. from Roskam
W_4_3 = 0.985; % Estimated climb to 40,000 ft W frac. + 0.005 for efficiency from Roskam and Dr. Barrett

% Breguet Range Equation for 3000 nmi range
R = 3000; % 3000 nmi range
M = 0.75; % Cruise speed as specified by objective function
c_j = 0.34; % TSFC (lbf/lbf hr) assuming better than current engines (0.35) but using greener fuels so efficiency is not super improved
L_D = 26; % L/D 26:1 assuming nlf P180 configuration
a_fts = 968.1; % Speed of sound (ft/s) at 40,000 ft
C_fts_kts = 3600/6076; % 3600 s/hr / 6076 ft/nmi
a = a_fts*C_fts_kts; % Convert speed of sound from ft/s to kts
V = M*a; % Velocity (kts)

W_5_4 = 1/exp((R*c_j)/(V*L_D)); % Rearranged B. Range Eqn. to find range W frac.

% W frac. continued
W_6_5 = 0.995; % Estimated descent W frac. +0.005 for efficiency from Roskam
W_7_6 = 1; % Estimated botched landing W frac. from Dr. Barrett
W_8_7 = 0.995; % Estimated small climb W frac. from Dr. Barrett

% Breguet Range Equation for divert to alternate 100 nmi away
R = 100; % 100 nmi range
M = 0.75; % Cruise speed as specified by objective function
a_fts = 1036.8; % 20,000 speed of sound ft/s (short climb altitude)
C_fts_kts = 3600/6076; % 3600 s/hr / 6076 ft/nmi
a = a_fts*C_fts_kts; % Convert speed of sound from ft/s to kts
V = M*a; % Divert to alternate velocity (kts)

W_9_8 = 1/exp((R*c_j)/(V*L_D)); % Rearranged B. Range Eqn. to find divert to alternate W frac.

% Breguet Endurance Equation for 45 min loiter after diversion per FAR 25
E_min = 45; % 45 min loiter endurance
E = E_min/60; % Convert to hr
W_10_9 = 1/exp((E*c_j)/(L_D)); % Estimated loiter W frac.

% W frac. continued
W_11_10 = 0.995; % Estimated descent W frac. +0.005 for efficiency from Roskam and Dr. Barrett
W_12_11 = 1; % Estimated landing W frac. from Dr. Barrett
W_13_12 = 0.995; % Estimated taxi and shutdown W frac. +0.003 for efficiency form Roskam and Dr. Barrett

% Mission Fuel Fraction Equation to find weight of fuel
M_ff = W_1_to*W_2_1*W_3_2*W_4_3*W_5_4*W_6_5*W_7_6*W_8_7*W_9_8*W_10_9*W_11_10*W_12_11*W_13_12;
W_f_used = (1-M_ff)*W_toguess1; % Weight (lbf) of fuel used

% Tentative Operating Empty Weight (lbf) Equation
W_oe_tent = W_toguess1 - W_f_used - W_crew

% Tentative Empty Weight (lbf) Equation and total fuel weight calculation
W_tfo = 0.01*W_f_used % Estimated trapped fuel and oil weight (lbf)
W_e_tent = W_oe_tent - W_tfo % Tentative empty weight (lbf)
W_f = W_f_used + W_tfo

