function [Eta_p, P_aero, T, Omega, alphaL] = PropDesign_SetPower(A,V_flt_kts, radius)

% Extract input parameters
Theta_0 = A(1);
Theta_1 = A(2);
Theta_tw = A(3);
Theta_r_twistExponent = A(4);
c_r = A(5);
P_desired = A(6);

% Define Parameters
nBlades = 3;
F_start = 1;
C_l_alpha = 0.108 / deg2rad(1); % lift curve slope NACA 0012 (1/radians) (Assumed Const)
C_D_0 = 0.01; % zero lift drag coefficient (Assumed Const)
d1 = 0.025; % (Assumed Const)
d2 = 0.65; % (Assumed Const)
aStall_minus_ao = deg2rad(17);
Omega_rpm = 1284;
rho_slug_ft3 = 0.002133;
% V_flt_kts = 75;
lambda = 0.8;
PowerError = 0.02;
nSegments = 100;
RootCuttout = 0.15;
% radius = 1; % m

% Convert units
Omega = Omega_rpm * 2 * pi / 60; % Convert RPM to rad/s
rho = rho_slug_ft3 * 515.379; % Convert slug/ft^3 to kg/m^3
V_flt = V_flt_kts * 0.514444; % Convert kts to m/s
dist_nonDim = 1 / nSegments;

% Precompute some constants
c_ave = (c_r + c_r * lambda) / 2;
sigma = nBlades * c_ave / (pi * radius);

% Required Unit Conversions
radius_eng = radius * 3.28084; % Convert meters to feet
c_r_eng = c_r * 3.28084; % Convert meters to feet
V_flt_eng = V_flt * 3.28084; % Convert m/s to ft/s
aStall_minus_ao_deg = 17; % Precomputed

% Preallocate arrays
alpha_limited = zeros(1, nSegments);
delta_C_t = zeros(1, nSegments);
delta_C_p_i = zeros(1, nSegments);
delta_C_p_o = zeros(1, nSegments);
alphaL = zeros(1, nSegments);


% Main iteration loop
i = 0;
while abs(PowerError) > 0.005
    i = i + 1;

    if PowerError > 0
        Omega = Omega - abs(Omega * PowerError)*.05;
    else 
        Omega = Omega + abs(Omega * PowerError)*.05;
    end

    % Initial Calculations
    c_ave_ft = c_ave * 3.28084; % Convert meters to feet
    V_tip = radius * Omega;
    J = V_flt / V_tip;

    % Begin Blade Segment Calculations
    for station = 1:nSegments
        [delta_C_t(station), delta_C_p_i(station), delta_C_p_o(station), alphaL(station)] = AnalyzeSegment(station, dist_nonDim, RootCuttout, nBlades, c_r, lambda, radius, V_flt, V_tip, Theta_0, Theta_1, Theta_tw, Theta_r_twistExponent, F_start, C_l_alpha, aStall_minus_ao, C_D_0, d1, d2);
    end

    % Combine Blade Segments
    C_T = sum(delta_C_t);
    C_P_i = sum(delta_C_p_i);
    C_P_o = sum(delta_C_p_o);
    C_P = C_P_i + C_P_o;

    % Post-processing
    T = C_T * (rho * pi * radius^2 * V_tip^2);
    P_shaft = C_P * (rho * pi * radius^2 * V_tip^3);
    P_aero = T * V_flt;
    Eta_p = P_aero / P_shaft;

    PowerError = (P_aero - P_desired) / P_desired;

    if i > 500
        disp('No convergence in 100 iterations');
        PowerError = 0; % end
    end
end

end

function [delta_C_t, delta_C_p_i, delta_C_p_o, alpha_limited] = AnalyzeSegment(station, dist_nonDim, RootCuttout, nBlades, c_r, lambda, radius, V_flt, V_tip, Theta_0, Theta_1, Theta_tw, Theta_r_twistExponent, F_start, C_l_alpha, aStall_minus_ao, C_D_0, d1, d2)

span = station * dist_nonDim;
midspan = span - dist_nonDim / 2;

solidity_local = (1 - heaviside(midspan - RootCuttout)) * 0.0001 + heaviside(midspan - RootCuttout) * nBlades * (c_r - c_r * span * (1 - lambda)) ./ (2 * pi * span * radius);

lambda_c = V_flt / V_tip;

Theta = Theta_0 + (Theta_1 + Theta_tw * midspan) .* midspan .^ Theta_r_twistExponent;

F = F_start;
lambda_r = ((solidity_local .* C_l_alpha ./ (16 .* F) - lambda_c ./ 2) .^ 2 + solidity_local .* C_l_alpha .* Theta .* midspan ./ 8) .^ 0.5 - (solidity_local .* C_l_alpha ./ 16 - lambda_c ./ 2);

for i = 2:5
    F_r = 2 ./ pi .* acos(exp(-nBlades ./ 2 .* (1 - midspan) ./ lambda_r));
    lambda_r = ((solidity_local .* C_l_alpha ./ (16 .* F_r) - lambda_c ./ 2) .^ 2 + solidity_local .* C_l_alpha .* Theta .* midspan ./ (8 .* F_r)) .^ 0.5 - (solidity_local .* C_l_alpha ./ (16 .* F_r) - lambda_c ./ 2);
end

alpha = Theta - atan2(lambda_r, midspan);
alpha_limited = min(alpha, aStall_minus_ao);

C_d = C_D_0 + d1 .* alpha_limited + d2 .* alpha_limited .^ 2;
C_l = C_l_alpha .* alpha_limited;

delta_C_t = solidity_local.* C_l_alpha / 2 .* (Theta .* midspan .^ 2 - lambda_r .* midspan) .* dist_nonDim;

% Handle special condition for the first station
if station == 1
    if alpha_limited * 180 / pi < 17
        delta_C_t = solidity_local * C_l_alpha / 2 * (Theta * midspan ^ 2 - lambda_r * midspan) * dist_nonDim;
    else
        delta_C_t = solidity_local * C_l_alpha / 2 * (aStall_minus_ao) * dist_nonDim;
    end
end

delta_C_p_i = lambda_r .* delta_C_t;
delta_C_p_o = solidity_local / 2 .* C_d .* midspan .^ 3 .* dist_nonDim;
    

end

