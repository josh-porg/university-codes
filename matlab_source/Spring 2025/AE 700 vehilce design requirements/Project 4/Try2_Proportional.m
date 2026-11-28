%% IR Seeker Missile Proportional Pursuit Simulation in XY Plane
% This simulation uses a proportional pursuit guidance law. The interceptor 
% (IR seeker) turns with a rate proportional to the error between its current 
% heading and the line-of-sight angle to the target. The turn rate is limited 
% to a maximum value.

clear; clc; close all;

%% Simulation Parameters
dt = 0.01;              % Time step [s]
T_total = 60;          % Total simulation time [s]
time = 0:dt:T_total;   % Time vector
numSteps = length(time);

%% Interceptor (Seeker) Parameters
S = [0, 0];                    % Initial interceptor position [m]
init_heading = pi/2;           % Initial heading: 90° (upward)
speed = 400;                   % Constant interceptor speed [m/s]
V = speed * [cos(init_heading), sin(init_heading)];  % Velocity vector

% Guidance parameters for proportional pursuit
Kp = .5;                      % Proportional gain for turn rate command [1/s]
max_turn_deg = 100;            % Maximum turn rate [deg/s]
max_turn_rate = deg2rad(max_turn_deg);  % Maximum turn rate [rad/s]

%% Target Parameters
% Target enters from the right and moves leftward.
x0_t = 300;             % Target initial x-position [m]
y0_t = 0;               % Target y-position (constant) [m]
Vt_mag = 200;           % Target speed [m/s] (target velocity vector is [-Vt_mag, 0])

%% IR Seeker Detector Parameters
Nd = 5;                 % Number of detectors (must be odd)
wd = 0.1;               % Total detector array width [m]
f = 0.05;               % Focal length [m]
delta_x_d = 0.01;       % Width of one detector element [m]
r_tgt_missile = 0.15;   % Radius of target missile [m] (30 cm diameter)
L_target = 0.6;         % Length of target missile [m] (e.g., 60 cm)

%% more guidance and fuzing parameters
swap_threshold = 2*Vt_mag/max_turn_rate;
collision_threshold = r_tgt_missile;

%% Fuze (Damage) Parameters
K = 90;                 % Number of fragments
phi = deg2rad(30);      % Fragmentation half-angle [rad] (example value)

%% Preallocate for Logging
S_history = zeros(numSteps,2);    % Interceptor path history
T_history = zeros(numSteps,2);    % Target path history
Pd_history = zeros(numSteps,1);     % Fuze probability over time
n_history = zeros(numSteps,1);      % Detector index output over time
Active_history = zeros(numSteps,1); % Illuminated detector count

%% Create Figure and Pause Button
hFig = figure;
hPlot = plot(NaN, NaN, 'b-', 'LineWidth', 2); % Interceptor path
hold on;
hTarget = plot(NaN, NaN, 'r-', 'LineWidth', 2);  % Target path
hInterceptor = plot(NaN, NaN, 'bo', 'MarkerFaceColor', 'b'); % Current interceptor
hTargetPos = plot(NaN, NaN, 'ro', 'MarkerFaceColor', 'r');    % Current target
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('IR Seeker Missile Proportional Pursuit Simulation');
legend('Interceptor Path', 'Target Path', 'Interceptor', 'Target');
grid on;

% Create a pause toggle button
setappdata(hFig, 'isPaused', false);
uicontrol('Style', 'togglebutton', ...
    'String', 'Pause', ...
    'Position', [20 20 60 30], ...
    'Callback', @(src, event)setappdata(hFig, 'isPaused', get(src,'Value')));

%% Simulation Loop
for i = 1:numSteps
    t = time(i);
    
    % Update Target Position: Target moves leftward
    %T_pos = [(-Vt_mag * t) + x0_t, y0_t];
    T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*sin(t/2) + .1*Vt_mag*sin(t*2) ) + y0_t];
    T_pos_est = [(-Vt_mag * (t+1)) + x0_t, y0_t]; % consider using this to correct for laggining the target
    
    % Call Detect_IR to obtain the detector index and active pixel count.
    % Function call: [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, 
    %           r_tgt_missile, L_target, Vt_mag, x0_t, y0_t)
    [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, r_tgt_missile, L_target, Vt_mag, x0_t, y0_t);
    
    % Compute line-of-sight (LOS) and proportional pursuit command:
    rel_pos = T_pos - S;                   % Relative position vector from interceptor to target
    desired_angle = atan2(rel_pos(2), rel_pos(1));  % LOS angle [rad]
    current_angle = atan2(V(2), V(1));       % Current interceptor heading [rad]
    
    % Compute the angle error between the desired LOS and current heading.
    angle_error = wrapToPi(desired_angle - current_angle);

    % adjust the persuite constant as you get closer top 1/dt as you approach the target
    Kp_adjusted = (1/dt) * 1/(1+exp(norm(S-T_pos)-swap_threshold)) + Kp * (1-1/(1+exp(norm(S-T_pos)-swap_threshold)))

    % Proportional pursuit guidance law: commanded turn rate = Kp * angle_error
    cmd_turn_rate = Kp_adjusted * angle_error;



    
    % % reduce the proportioanlity constant as a function of distance to
    % target
    % if norm(S-T_pos) < swap_threshold        
    %     cmd_turn_rate = angle_error/dt
    % end
    

    % Enforce the maximum turn rate limit.
    if abs(cmd_turn_rate) > max_turn_rate
        cmd_turn_rate = sign(cmd_turn_rate) * max_turn_rate;
    end
    
    % Update interceptor heading using the commanded turn rate.
    new_angle = current_angle + cmd_turn_rate * dt;
    V = speed * [cos(new_angle), sin(new_angle)];
    
    % Update interceptor position using new velocity.
    S = S + V * dt;
    
    % Calculate fuze/damage probability using the Fuze function:
    % Function call: [Pd] = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag)
    Pd = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag)
    
    % Log data for plotting and review
    S_history(i,:) = S;
    T_history(i,:) = T_pos;
    Pd_history(i) = Pd;
    n_history(i) = n;
    Active_history(i) = Active;
    
    % Update plot for interceptor and target paths and current positions
    set(hPlot, 'XData', S_history(1:i,1), 'YData', S_history(1:i,2));
    set(hTarget, 'XData', T_history(1:i,1), 'YData', T_history(1:i,2));
    set(hInterceptor, 'XData', S(1), 'YData', S(2));
    set(hTargetPos, 'XData', T_pos(1), 'YData', T_pos(2));
    drawnow;
    % pause(0.1)

    
    

    % Check for collision: if the distance is less than the threshold, stop simulation.
    if norm(S - T_pos) < collision_threshold
        collisionDetected = true;
        disp('Contact detected! Exiting simulation.');
        break;  % Exit simulation loop.
    end

    % Check for pause toggle
    while getappdata(hFig, 'isPaused')
        pause(0.1);
    end
end

%% Final Output Display
disp(['Final Probability of Damage (Pd): ' num2str(Pd)]);

% Plot Probability of Damage over time in a separate figure.
figure;
plot(time, Pd_history, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Probability of Damage (Pd)');
title('Fuze Function Output Over Time');
grid on;

%% --- Local Functions ---
% The following functions are defined as subfunctions appended at the end of this file.

function [n,Active] = Detect_IR(V,S,Nd,wd,f,delta_x_d,time_elapsed,...
    r_tgt_missile,L,Vt_mag,x0_t,y0_t)
% This function calculates the detector index (n) and the number of illuminated
% detectors (Active) for the IR seeker using the current interceptor and target
% parameters.

% Compute target's current position relative to some fixed frame.
St = [(-Vt_mag * time_elapsed) + x0_t, y0_t]; 
% Relative target position with respect to the interceptor.
St_prime = St - S;

% Compute the angle between the interceptor's velocity V and the line from
% interceptor to target.
theta_need = acos(dot(V,St_prime)/(norm(V)*norm(St_prime)));
% Determine sign for the needed turn angle using the cross product.
check = cross([St_prime 0], [V 0]); 
if check(3) > 0
    theta_need = -theta_need;
end

% Compute maximum allowable angle based on the detector array geometry.
theta_max = atan(wd/(2*f));

% If the required turn angle exceeds the detector's field of view, no detection.
if abs(theta_need) >= theta_max 
    n = 0;
    Active = 0;
    fprintf('No Detect \n');
    return
else
    % Map theta_need to the detector element index.
    n = (theta_need/theta_max) * (Nd/2);
    if theta_need <= 0
        n = floor(n) + 1; 
    else
        n = floor(n);
    end
end

% Calculate the projected area onto the detector.  
Vt = [-Vt_mag 0]; % Target velocity vector.
top = cross([V 0], [Vt 0]);
sin_theta_off = top(3)/(norm(V)*norm(Vt));
cos_psi = sin_theta_off;
A = 2 * r_tgt_missile * L;
Ae = A * cos_psi;
Le = L * cos_psi;
if Ae < pi * r_tgt_missile^2
    Le = 2 * r_tgt_missile;
end

% Determine the illumination scale and count the number of active detector pixels.
x_illum = (Le * f) / norm(St_prime);
Active = 2 * (floor(x_illum/delta_x_d) + 1) - 1;
end

function [Pd] = Fuze(V, r_tgt_missile, L, K, S, phi, time_elapsed, x0_t, y0_t, Vt_mag)
% This function calculates the probability of damage (Pd) based on the current 
% interceptor state and target parameters using an exponential model.

% Compute target's current position.
St = [(-Vt_mag * time_elapsed) + x0_t, y0_t]; 
% Relative target position with respect to the interceptor.
St_prime = St - S;
Range = norm(St_prime);

% Compute the off-axis angle from the interceptor velocity to target velocity.
Vt = [-Vt_mag 0]; 
top = cross([V 0], [Vt 0]);
sin_theta_off = top(3)/(norm(V)*norm(Vt));
cos_psi = sin_theta_off;
A = 2 * r_tgt_missile * L; 
Ae = A * cos_psi;
if Ae < pi * r_tgt_missile^2
    Ae = pi * r_tgt_missile^2;
end

% Compute the effective spreading area at the target range.
As = 2 * pi * (Range^2) * (1 - cos(phi));
% Exponential fuze/damage probability equation.
Pd = 1 - exp(-K * (Ae/As));
end
