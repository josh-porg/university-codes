%% IR Seeker Missile Pure Pursuit with Collision and Explosion Simulation in XY Plane
% This simulation uses pure pursuit guidance. The interceptor always steers
% directly toward the target’s current location. When the missile makes contact
% with the target (i.e. when the center-to-center distance is below a threshold),
% the simulation stops and an explosion effect is animated.
%
% The functions Detect_IR and Fuze are appended at the end of this file.

clear; clc; close all;

%% Simulation Parameters
dt = 0.01;              % Time step [s]
T_total = 60;          % Total simulation time [s]
time = 0:dt:T_total;   % Time vector
numSteps = length(time);

%% Interceptor (Seeker) Parameters
S = [0, 0];            % Initial interceptor position [m]
init_heading = pi/2;   % Initial heading (90° upward) in radians
speed = 400;           % Interceptor speed [m/s]
V = speed * [cos(init_heading), sin(init_heading)];  % Initial velocity vector

max_turn_deg = 100;                % Maximum turn rate [deg/s]
max_turn_rate = deg2rad(max_turn_deg);  % Maximum turn rate in rad/s
max_angle_change = max_turn_rate * dt;    % Maximum heading change per time step [rad]

%% Target Parameters
% Target enters from the right and moves leftward.
x0_t = 300;            % Target initial x-position [m]
y0_t = 0;              % Target y-position [m]
Vt_mag = 200;          % Target speed [m/s] (target moves leftward)
target_amplitude = 50;  % Vertical oscillation amplitude (m)
target_period = 20;     % Oscillation period (s)

%% IR Seeker Detector Parameters
Nd = 5;                % Number of detectors (must be odd)
wd = 0.1;              % Total width of the detector array [m]
f = 0.05;              % Focal length [m]
delta_x_d = 0.01;      % Width of one detector element [m]
r_tgt_missile = 0.15;  % Target missile radius [m] (30 cm diameter)
L_target = 0.6;        % Target missile length [m]

%% Fuze (Damage) Parameters
K = 90;                % Number of fragments
phi = deg2rad(30);     % Fragmentation half-angle [rad]

%% Collision Parameters
% We stop the simulation when the center-to-center distance is below this threshold.
collision_threshold = r_tgt_missile;  

%% Preallocate Logging Arrays
S_history = zeros(numSteps,2);    % Interceptor path history
T_history = zeros(numSteps,2);      % Target path history
Pd_history = zeros(numSteps,1);     % Fuze probability history
n_history = zeros(numSteps,1);      % Detector index history
Active_history = zeros(numSteps,1); % Illuminated detector count history

%% Create Figure and Pause Button
hFig = figure;
hPlot = plot(NaN, NaN, 'b-', 'LineWidth', 2);  % Interceptor path
hold on;
hTarget = plot(NaN, NaN, 'r-', 'LineWidth', 2);  % Target path
hInterceptor = plot(NaN, NaN, 'bo', 'MarkerFaceColor', 'b');  % Current interceptor position
hTargetPos = plot(NaN, NaN, 'ro', 'MarkerFaceColor', 'r');     % Current target position
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('IR Seeker Missile Pure Pursuit Simulation');
legend('Interceptor Path','Target Path','Interceptor','Target','Location','Best');
grid on;

% Create a toggle button for pausing the simulation
setappdata(hFig, 'isPaused', false);
uicontrol('Style', 'togglebutton', ...
    'String', 'Pause', ...
    'Position', [20 20 60 30], ...
    'Callback', @(src, event)setappdata(hFig, 'isPaused', get(src,'Value')));

%% Simulation Loop
collisionDetected = false;  % Flag to exit the simulation loop upon impact
for i = 1:numSteps
    t = time(i);
    
    % Update target position (target moves leftward).
    T_pos = [(-Vt_mag * t) + x0_t, target_amplitude * sin(2*pi*t/target_period)];
    T_pos = [(-Vt_mag*.5*cos(t) * t) + x0_t, (-Vt_mag*.5*sin(t) * t) + y0_t];
    
    % Call Detect_IR to obtain the detector index and active pixel count.
    % Function call:
    %   [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, ...
    %                             r_tgt_missile, L_target, Vt_mag, x0_t, y0_t)
    try
        [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, r_tgt_missile, L_target, Vt_mag, x0_t, target_amplitude, target_period);
    catch
    end
    % PURE PURSUIT GUIDANCE:
    % Compute the line-of-sight (LOS) angle.
    rel_pos = T_pos - S;            % Relative vector from interceptor to target
    desired_angle = atan2(rel_pos(2), rel_pos(1));  % LOS angle [rad]
    current_angle = atan2(V(2), V(1));             % Current interceptor heading [rad]
    
    % Compute angle error and update heading.
    angle_error = wrapToPi(desired_angle - current_angle);
    if abs(angle_error) <= max_angle_change
        new_angle = desired_angle;
    else
        new_angle = current_angle + sign(angle_error)*max_angle_change;
    end
    V = speed * [cos(new_angle), sin(new_angle)]; % Update velocity vector
    
    % Update interceptor position.
    S = S + V * dt;
    
    % Calculate fuze/damage probability.
    % Function call:
    %   [Pd] = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag)
    Pd = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, target_amplitude, target_period, Vt_mag);
    
    % Log simulation data.
    S_history(i,:) = S;
    T_history(i,:) = T_pos;
    Pd_history(i) = Pd;
    n_history(i) = n;
    Active_history(i) = Active;
    
    % Update plots for interceptor and target.
    set(hPlot, 'XData', S_history(1:i,1), 'YData', S_history(1:i,2));
    set(hTarget, 'XData', T_history(1:i,1), 'YData', T_history(1:i,2));
    set(hInterceptor, 'XData', S(1), 'YData', S(2));
    set(hTargetPos, 'XData', T_pos(1), 'YData', T_pos(2));
    drawnow;
    
    % Check for collision: if the distance is less than the threshold, stop simulation.
    if norm(S - T_pos) < collision_threshold
        collisionDetected = true;
        disp('Contact detected! Exiting simulation.');
        break;  % Exit simulation loop.
    end
    
    % Pause handling.
    while getappdata(hFig, 'isPaused')
        pause(0.1);
    end
end

% End of Simulation
if ~collisionDetected
    disp('60 seconds reached without contact. Simulation ended.');
end

%% If a collision was detected, animate an explosion.
if collisionDetected
    explosion_center = T_pos;  % Use the target's (or interceptor's) position.
    explosion_final_radius = 10;  % Final explosion radius [m] (adjust for visibility)
    num_frames = 15;
    theta = linspace(0,2*pi,100);
    for r = linspace(0, explosion_final_radius, num_frames)
         x_ex = explosion_center(1) + r*cos(theta);
         y_ex = explosion_center(2) + r*sin(theta);
         % Draw explosion as a filled red circle with no edge.
         h_explosion = fill(x_ex, y_ex, 'r','EdgeColor','none');
         drawnow;
         pause(0.05);
         if r < explosion_final_radius
             delete(h_explosion);
         end
    end
    % Leave the final explosion visible and add a text label.
    text(explosion_center(1), explosion_center(2), 'BOOM!', 'HorizontalAlignment','center', ...
         'VerticalAlignment','middle', 'FontSize',14, 'Color','w', 'FontWeight','bold');
end

%% Final Output Display
disp(['Final Probability of Damage (Pd): ' num2str(Pd)]);

% Plot the fuze probability over time in a separate figure.
figure;
plot(time(1:i), Pd_history(1:i), 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Probability of Damage (Pd)');
title('Fuze Function Output Over Time');
grid on;

%% --- Local Functions ---

function [n,Active] = Detect_IR(V,S,Nd,wd,f,delta_x_d,time_elapsed,...
    r_tgt_missile,L,Vt_mag,x0_t,target_amplitude,target_period)
% Detect_IR calculates the detector index (n) and number of illuminated
% detector pixels (Active) for the IR seeker based on the interceptor and target state.
%
% Inputs:
%   V            - Interceptor velocity vector [m/s]
%   S            - Interceptor position vector [m]
%   Nd           - Number of detectors (odd)
%   wd           - Total width of detector array [m]
%   f            - Focal length [m]
%   delta_x_d  - Width of one detector element [m]
%   time_elapsed - Current simulation time [s]
%   r_tgt_missile - Target missile radius [m]
%   L            - Target missile length [m]
%   Vt_mag       - Target speed [m/s]
%   x0_t, y0_t   - Target initial positions [m]
%
% Outputs:
%   n      - Detector index (relative to center; 0 indicates no detection)
%   Active - Number of illuminated detector pixels

% Compute target's current position.
St = [(-Vt_mag*time_elapsed) + x0_t, target_amplitude*sin(2*pi*time_elapsed/target_period)];
St_prime = St - S;   % Relative target position

% Compute angle between interceptor's velocity and the line to target.
theta_need = acos(dot(V,St_prime)/(norm(V)*norm(St_prime)));
check = cross([St_prime 0],[V 0]);
if check(3) > 0
    theta_need = -theta_need;
end

% Maximum allowable deviation based on detector geometry.
theta_max = atan(wd/(2*f));

if abs(theta_need) >= theta_max
    n = 0;
    Active = 0;
    fprintf('No Detect \n');
    return
else
    n = (theta_need/theta_max) * (Nd/2);
    if theta_need <= 0
        n = floor(n) + 1;
    else
        n = floor(n);
    end
end

% Calculate projected area to determine active pixels.
Vt = [-Vt_mag 0];
top = cross([V 0],[Vt 0]);
sin_theta_off = top(3)/(norm(V)*norm(Vt));
cos_psi = sin_theta_off;
A = 2 * r_tgt_missile * L;
Ae = A * cos_psi;
Le = L * cos_psi;
if Ae < pi * r_tgt_missile^2
    Le = 2 * r_tgt_missile;
end

x_illum = (Le * f) / norm(St_prime);
Active = 2 * (floor(x_illum/delta_x_d) + 1) - 1;
end

function [Pd] = Fuze(V, r_tgt_missile, L, K, S, phi, time_elapsed, x0_t, target_amplitude, target_period, Vt_mag)
% Fuze calculates the probability of damage (Pd) based on the interceptor state
% and target parameters using an exponential model.
%
% Inputs:
%   V             - Interceptor velocity vector [m/s]
%   r_tgt_missile - Target missile radius [m]
%   L             - Target missile length [m]
%   K             - Number of fragments
%   S             - Interceptor position vector [m]
%   phi           - Fragmentation half-angle [rad]
%   time_elapsed  - Current simulation time [s]
%   x0_t, y0_t    - Target initial positions [m]
%   Vt_mag        - Target speed [m/s]
%
% Output:
%   Pd - Probability of damage

% Compute target's current position on a curved trajectory.
St = [(-Vt_mag*time_elapsed) + x0_t, target_amplitude*sin(2*pi*time_elapsed/target_period)];
St_prime = St - S;
Range = norm(St_prime);

Vt = [-Vt_mag 0];
top = cross([V 0],[Vt 0]);
sin_theta_off = top(3)/(norm(V)*norm(Vt));
cos_psi = sin_theta_off;
A = 2 * r_tgt_missile * L;
Ae = A * cos_psi;
if Ae < pi * r_tgt_missile^2
    Ae = pi * r_tgt_missile^2;
end

As = 2 * pi * (Range^2) * (1 - cos(phi));
Pd = 1 - exp(-K * (Ae/As));
end
