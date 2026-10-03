%% IR Seeker Missile Proportional Pursuit Simulation in XY Plane
% This simulation uses a proportional pursuit guidance law. The interceptor 
% (IR seeker) turns with a rate proportional to the error between its current 
% heading and the line-of-sight angle to the target. The turn rate is limited 
% to a maximum value.

clear; clc; close all;

%% Simulation Parameters
dt = 0.001;              % Time step [s]
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
Nd = 2.1998e+04;        % Number of detectors (must be odd)
wd = 0.2750;               % Total detector array width [m]
f = 0.0030;               % Focal length [m]
delta_x_d = 1.2500e-05;       % Width of one detector element [m]
r_tgt_missile = 0.15;   % Radius of target missile [m] (30 cm diameter)
L_target = 0.6;         % Length of target missile [m] (e.g., 60 cm)

%% IR Seeker Detector Parameters
Nd = 7.3327e+07;        % Number of detectors (must be odd)
wd = 0.2750;               % Total detector array width [m]
f = 0.0030;               % Focal length [m]
delta_x_d = 3.7500e-08;       % Width of one detector element [m]
r_tgt_missile = 0.15;   % Radius of target missile [m] (30 cm diameter)
L_target = 0.6;         % Length of target missile [m] (e.g., 60 cm)

%% more guidance and fuzing parameters
swap_threshold = 2*Vt_mag/max_turn_rate *100;
swap_threshold = wd/delta_x_d * (.00001);
collision_threshold = r_tgt_missile;

%% porpotional guidance variable initialization
lambda_dot_tgt = 0; % target angular rate [rad/s]
V_c = 0; % target closure velocity [rad/s]
r_to_tgt_old = 1000; % target previous radial distance [m]
lambda_tgt_old = -pi/4; % target previous angular position [rad]
T_pos_est_old = [0,0]
Kp_adjusted_old = .5

%% Fuze (Damage) Parameters
K = 90;                 % Number of fragments
phi = deg2rad(30);      % Fragmentation half-angle [rad] (example value)

%% kalman filter parameters
% Initialize state variables for the moving target
x_target = 50;        % Initial x position of the target
y_target = 50;        % Initial y position of the target
vx_target = 5;        % Initial velocity of the target in x
vy_target = 2;        % Initial velocity of the target in y

x_missile = 0;       % Initial x position of the missile
y_missile = 0;       % Initial y position of the missile
vx_missile = 10;     % Initial velocity of the missile in xH
vy_missile = 5;      % Initial velocity of the missile in y
% Process and measurement noise
Q = eye(4);          % Process noise covariance
R = 10*eye(2);       % Measurement noise covariance

% Measurement noise covariance
R = 10 * eye(2);       % Measurement noise covariance for target position
R_IR = 1 * eye(2);      % Measurement noise covariance for IR sensor (adjust as needed)

% Initialize Kalman filter parameters
A = eye(4);          % State transition matrix
B = eye(4);          % Control matrix
H = eye(2,4);        % Measurement matrix

% Initial state estimate
x_hat = [x_missile; y_missile; vx_missile; vy_missile];

% Initialize state covariance matrix
P = eye(4);

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
hTargetPos_est = plot(NaN, NaN, 'gx', 'MarkerFaceColor', 'g');    % Current target estimate
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


%stochastic evasion
% Parameters
n = length(time); % Number of points
cutoffFreq = .001; % Cutoff frequency for low-pass filter
% Generate random noise
noise = randn(1, n);

% Design low-pass filter
[b, a] = butter(2, cutoffFreq, 'low');

% Apply low-pass filter
smoothNoise = filtfilt(b, a, noise);

%% Simulation Loop
for i = 1:numSteps
    t = time(i);

    
    
    
    
    % Update Target Position: Target moves leftward
    T_pos = [(-Vt_mag * t) + x0_t, y0_t];
    T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*sin(t/2) + .1*Vt_mag*sin(t*2) ) + y0_t];
    T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t];
    %T_pos_est = [(-Vt_mag * (t+1)) + x0_t, y0_t]; % consider using this to correct for laggining the target
    
    % Call Detect_IR to obtain the detector index and active pixel count.
    % Function call: [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, 
    %           r_tgt_missile, L_target, Vt_mag, x0_t, y0_t)
    [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, r_tgt_missile, L_target, Vt_mag, x0_t, y0_t);

    % must now estimate position and velocity
    % let's try something really simple.
    % assume the target is the average of the minimum and maximum exposed
    % lengths
    A_tgt = ((L_target * r_tgt_missile) + (pi * r_tgt_missile^2)) /2;
    l_tgt_apparent = (L_target + r_tgt_missile) /2; % TODO: use angular size (estimate from its angle)

    % now we find the distance from the sensor using its size
    r_to_tgt = l_tgt_apparent / (Active * delta_x_d); % from similar triangles

    % now we want its angular position
    lambda_tgt = wrapTo2Pi(-atan2(n * delta_x_d, f));

    %compute derivatives with central difference
    lambda_dot = (lambda_tgt - lambda_tgt_old) /dt; % agular rate
    V_c = (r_to_tgt - r_to_tgt_old) / dt; % closure velocity

    % if isinf(r_to_tgt)
    %     r_to_tgt = r_to_tgt_old + V_c;
    %     lambda_tgt = lambda_tgt_old + lambda_dot_tgt;
    % end

    % if isinf(r_to_tgt)
    %     T_pos_est = T_pos_est_old
    % end

    if (r_to_tgt - r_to_tgt_old) / ((r_to_tgt_old)) > 1
        r_to_tgt = (r_to_tgt+r_to_tgt_old)/2
    end

     % estimate target position
    if(isnan(lambda_tgt) || isnan(r_to_tgt) || isinf(r_to_tgt))
        current_angle = atan2(V(2), V(1))
        T_pos_est = S + [r_to_tgt_old *cos(lambda_tgt_old + current_angle), r_to_tgt_old *sin(lambda_tgt_old + current_angle)]
    else
        % estimate target position   
        T_pos_est = S + [r_to_tgt *cos(lambda_tgt), r_to_tgt *sin(lambda_tgt)];

        % % Kalman filter prediction
        % x_hat_minus = A * x_hat;
        % P_minus = A * P * A' + Q;
        % 
        % % Kalman filter update
        % K_kf = P_minus * H' / (H * P_minus * H' + R_IR);
        % x_hat = x_hat_minus + K_kf * (T_pos_est - H * x_hat_minus);
        % P = (eye(4) - K_kf * H) * P_minus;
        % 
        % T_pos_est = x_hat(1:2)
    end

    

    if isnan(T_pos_est)
        T_pos_est = S + [r_to_tgt_old *cos(lambda_tgt_old), r_to_tgt_old *sin(lambda_tgt_old)];
    elseif(~isnan(lambda_tgt) && ~isnan(r_to_tgt) && ~isinf(r_to_tgt))
        lambda_tgt_old = lambda_tgt;
        r_to_tgt_old = r_to_tgt;
    end

    T_pos_est_old = T_pos_est;
    

    % % Kalman filter prediction
    % x_hat_minus = A * x_hat;
    % P_minus = A * P * A' + Q;
    % 
    % % Kalman filter update
    % K_kf = P_minus * H' / (H * P_minus * H' + R_IR);
    % x_hat = x_hat_minus + K_kf * (T_pos_est - H * x_hat_minus);
    % P = (eye(4) - K_kf * H) * P_minus;
    % 
    % T_pos_est = x_hat(1:2)

    



    
    % Compute line-of-sight (LOS) and proportional pursuit command:
    rel_pos = T_pos - S;                   % Relative position vector from interceptor to target
    rel_pos = T_pos_est(1:2) - S; % use estimated position for guidance
    desired_angle = atan2(rel_pos(2), rel_pos(1));  % LOS angle [rad]
    current_angle = atan2(V(2), V(1));       % Current interceptor heading [rad]
    
    % Compute the angle error between the desired LOS and current heading.
    angle_error = wrapToPi(desired_angle - current_angle);

    % adjust the persuite constant as you get closer top 1/dt as you approach the target
    Kp_adjusted = (1/dt) * 1/(1+exp(norm(S-T_pos_est)-swap_threshold)) + Kp * (1-1/(1+exp(norm(S-T_pos_est)-swap_threshold)));
    Kp_adjusted = Kp * 1/(1+exp(Active-swap_threshold)) + (1/dt) * (1-1/(1+exp(Active-swap_threshold)))

    
    Kp_adjusted = (Kp_adjusted*1/4 + Kp_adjusted_old*3/4);

    Kp_adjusted_old = Kp_adjusted

    % Proportional pursuit guidance law: commanded turn rate = Kp * angle_error
    cmd_turn_rate = Kp_adjusted * angle_error;

    % % reduce the proportioanlity constant as a function of distance to
    % target
    % if norm(S-T_pos) < swap_threshold        
    %     cmd_turn_rate = angle_error/dt
    % end

    if(r_to_tgt == Inf)
        cmd_turn_rate = pi*50*(lambda_tgt_old);
        cmd_turn_rate = pi*50* sign(lambda_tgt_old);
    end
    

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
    Pd = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag, T_pos)
    
    % Log data for plotting and review
    S_history(i,:) = S;
    T_history(i,:) = T_pos;
    Pd_history(i) = Pd;
    n_history(i) = n;
    Active_history(i) = Active;
    
    if mod(i,50) == 0
        % Update plot for interceptor and target paths and current positions
        set(hPlot, 'XData', S_history(1:i,1), 'YData', S_history(1:i,2));
        set(hTarget, 'XData', T_history(1:i,1), 'YData', T_history(1:i,2));
        set(hInterceptor, 'XData', S(1), 'YData', S(2));
        set(hTargetPos, 'XData', T_pos(1), 'YData', T_pos(2));
        %set(hTargetPos_est, 'XData', T_pos_est(1), 'YData', T_pos_est(2))
        drawnow;
        % pause(0.1)
    end

    
    

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

function [Pd] = Fuze(V, r_tgt_missile, L, K, S, phi, time_elapsed, x0_t, y0_t, Vt_mag, T_pos)
% This function calculates the probability of damage (Pd) based on the current 
% interceptor state and target parameters using an exponential model.

% Compute target's current position.
St = [(-Vt_mag * time_elapsed) + x0_t, y0_t];
St = [(-Vt_mag*time_elapsed ) + x0_t, (Vt_mag*sin(time_elapsed/2) + .1*Vt_mag*sin(time_elapsed*2) ) + y0_t];
St = T_pos;
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


function noise = rands(sz, sig)
    if nargin < 2
        sig = 0.1; % Default Gaussian width
    end
    % Generate random numbers
    randomNumbers = randn(sz);
    % Smooth the random numbers with a Gaussian kernel
    noise = conv2(randomNumbers, fspecial('gaussian', [5 5], sig), 'same');
end