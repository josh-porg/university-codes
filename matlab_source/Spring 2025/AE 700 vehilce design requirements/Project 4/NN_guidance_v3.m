% -------------------------------------------------------------------------
% 1. Create a feedforward network with two hidden layers of 32 neurons each
% -------------------------------------------------------------------------
hiddenSizes = [32, 32];  
net = feedforwardnet(hiddenSizes);         
net = configure(net, rand(2,1), 0);
net.trainParam.showWindow = false;

% 2. Turn off automatic input/output normalization
netTemplate.inputs{1}.processFcns  = {};               % remove removeconstantrows & mapminmax
netTemplate.outputs{2}.processFcns = {};               % same for output

% By default, feedforwardnet expects inputs of size [numInputs x numSamples].
% Our inputs are [n; Active] from Detect_IR, so we configure the network
% without training data purely to set input/output dimensions:
dummyInput  = rand(2,1);                    % 2×1 dummy sample
dummyTarget = zeros(1,1);                   % 1×1 dummy target
net = configure(net, dummyInput, dummyTarget); 

% Choose a training function (if you later train supervised):
net.trainFcn = 'trainlm';    % Levenberg–Marquardt backprop (default)
net.performFcn = 'mse';      % Mean‑squared error
% Optionally turn off training GUI:
net.trainParam.showWindow = false;

% 3. Determine decision‐vector length
wb0   = getwb(net);           % initial weights/biases
nVars = numel(wb0);                   % total number of GA vars


%% attempt two custom NN
% 1) Build initial w0
layerSizes = [32,32];
nVars = layerSizes(1)*2 + layerSizes(1) ...     % weights1 + bias1
      + layerSizes(2)*layerSizes(1) + layerSizes(2) ... % weights2 + b2
      + layerSizes(2) + 1;                       % W3 + b3

w0 = randn(1, nVars) * 0.1;    % small random init


% 4. Set GA bounds (e.g. ±1 around zero)
lb = -1*ones(1, nVars);
ub =  1*ones(1, nVars);

% 5. Create a fitness-function handle (passing extra args)
fitnessFcn = @(x) computeFitness(x, net);

% 6. Configure GA options
popSize = 1000; % how many individual per generation
elitePercent = .10; % save this fraction of population (the best)

opts = optimoptions('ga', ...
    'PopulationSize', popSize, ...
    'MaxGenerations', 70, ...
    'UseParallel', true, ...
    'EliteCount', ceil(elitePercent*popSize), ...
    'Display', 'iter');

% 7. Run GA!
[xBest, fvalBest] = ga(fitnessFcn, nVars, [], [], [], [], lb, ub, [], opts);

% 8. Build the optimized network

bestNN = @(n,Active) mySimpleNN(xBest, [n;Active], layerSizes);
% netOptimized = setwb(netTemplate, xBest);
% disp(['Optimized max Pd = ', num2str(1 - fvalBest)]);

maxPd = simulateAndPlot(xBest, net);


function Pd_avg = computeFitness(x, net)
    Pd_avg = 0;
    n_inters = 10;
    % do it 3 times because the running is random
    for itteration = 1:n_inters
        Pd_avg = Pd_avg + simulate(x, net);
    end

    Pd_avg = Pd_avg/n_inters;
    close all;
end


function maxPd = simulate(x, net)
    %% assemble the net
    % Assign weights/biases to a fresh copy of the NN
    nnFcn = @(n,Active) mySimpleNN(x, [n;Active], [3,3]);
    %net = setwb(net, x);

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
    Nd = 7.3327e+06;        % Number of detectors (must be odd)
    wd = 0.2750;               % Total detector array width [m]
    f = 0.0030;               % Focal length [m]
    delta_x_d = 3.7500e-08;       % Width of one detector element [m]
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
    % hFig = figure;
    % hPlot = plot(NaN, NaN, 'b-', 'LineWidth', 2); % Interceptor path
    % hold on;
    % hTarget = plot(NaN, NaN, 'r-', 'LineWidth', 2);  % Target path
    % hInterceptor = plot(NaN, NaN, 'bo', 'MarkerFaceColor', 'b'); % Current interceptor
    % hTargetPos = plot(NaN, NaN, 'ro', 'MarkerFaceColor', 'r');    % Current target
    % xlabel('X Position (m)');
    % ylabel('Y Position (m)');
    % title('IR Seeker Missile Proportional Pursuit Simulation');
    % legend('Interceptor Path', 'Target Path', 'Interceptor', 'Target');
    % grid on;
    
    % Create a pause toggle button
    % setappdata(hFig, 'isPaused', false);
    % uicontrol('Style', 'togglebutton', ...
    %     'String', 'Pause', ...
    %     'Position', [20 20 60 30], ...
    %     'Callback', @(src, event)setappdata(hFig, 'isPaused', get(src,'Value')));

    % make missile run in random direction
    run_angle = rand()*2*pi;

    Rot = [cos(run_angle) -sin(run_angle); 
     sin(run_angle)  cos(run_angle)];

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
        %T_pos = [(-Vt_mag * t) + x0_t, y0_t];
        T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*sin(t/2) + .1*Vt_mag*sin(t*2) ) + y0_t];
        T_pos = [((-Vt_mag*t ) + x0_t )*cos(run_angle) - (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t*sin(run_angle), ((-Vt_mag*t ) + x0_t )*cos(run_angle) + (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t *sin(run_angle)];
        T_pos = (Rot * T_pos')';
        %T_pos_est = [(-Vt_mag * (t+1)) + x0_t, y0_t]; % consider using this to correct for laggining the target
        
        % Call Detect_IR to obtain the detector index and active pixel count.
        % Function call: [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, 
        %           r_tgt_missile, L_target, Vt_mag, x0_t, y0_t)
        [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, r_tgt_missile, L_target, Vt_mag, x0_t, y0_t);
    
        % --- 2. NN turn‑rate command
        inp = [n; Active];                   % 2×1
        %cmd_turn_rate = net(inp);                 % network output (1×1)
        cmd_turn_rate = nnFcn(n, Active);  % where nnFcn = @(n,Active) mySimpleNN(x, [n;Active], layerSizes)


        % --- 3. Limit and apply turn rate
        % Enforce the maximum turn rate limit.
        if abs(cmd_turn_rate) > max_turn_rate
            cmd_turn_rate = sign(cmd_turn_rate) * max_turn_rate;
        end
        
        % Update interceptor heading using the commanded turn rate.
        current_angle = atan2(V(2), V(1));       % Current interceptor heading [rad]
        new_angle = current_angle + cmd_turn_rate * dt;
        V = speed * [cos(new_angle), sin(new_angle)];
        
        % Update interceptor position using new velocity.
        S = S + V * dt;
        
        % Calculate fuze/damage probability using the Fuze function:
        % Function call: [Pd] = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag)
        Pd = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag, T_pos);
        
        % Log data for plotting and review
        S_history(i,:) = S;
        T_history(i,:) = T_pos;
        Pd_history(i) = Pd;
        n_history(i) = n;
        Active_history(i) = Active;
        
        % if mod(i,50) == 0
        %     % Update plot for interceptor and target paths and current positions
        %     set(hPlot, 'XData', S_history(1:i,1), 'YData', S_history(1:i,2));
        %     set(hTarget, 'XData', T_history(1:i,1), 'YData', T_history(1:i,2));
        %     set(hInterceptor, 'XData', S(1), 'YData', S(2));
        %     set(hTargetPos, 'XData', T_pos(1), 'YData', T_pos(2));
        %     drawnow;
        %     % pause(0.1)
        % end
    
        
        % Early exit on intercept (optional)
        if Pd > 0.9
            break;
        end
    
        % Check for collision: if the distance is less than the threshold, stop simulation.
        if norm(S - T_pos) < collision_threshold
            collisionDetected = true;
            disp('Contact detected! Exiting simulation.');
            break;  % Exit simulation loop.
        end
        
    
        % % Check for pause toggle
        % while getappdata(hFig, 'isPaused')
        %     pause(0.1);
        % end
    end

    maxPd = max(Pd_history);

    maxPd = -maxPd; % negative is good for minimization
    close
end

%% mySimpleNN.m
% A feed‑forward network with two hidden layers you control by a weight vector.
% Usage: y = mySimpleNN(w, inp, layerSizes)
%
%   w          : 1×N vector of all weights & biases
%   inp        : 2×1 input vector [n; Active]
%   layerSizes : [H1, H2] sizes of hidden layers
%   y          : scalar network output

function y = mySimpleNN(w, inp, layerSizes)
    % Unpack layer sizes
    H1 = layerSizes(1);
    H2 = layerSizes(2);

    % Initialize index
    idx = 0;

    % Extract weights1 and reshape to H1 x 2
    nweights1 = H1 * 2;
    weights1 = reshape(w(idx + 1:idx + nweights1), H1, 2);
    idx = idx + nweights1;

    % Extract bias1 and reshape to H1 x 1
    bias1 = reshape(w(idx + 1:idx + H1), H1, 1);
    idx = idx + H1;

    % Extract weights2 and reshape to H2 x H1
    nweights2 = H2 * H1;
    weights2 = reshape(w(idx + 1:idx + nweights2), H2, H1);
    idx = idx + nweights2;

    % Extract b2 and reshape to H2 x 1
    b2 = reshape(w(idx + 1:idx + H2), H2, 1);
    idx = idx + H2;

    % Extract W3 and reshape to 1 x H2
    W3 = reshape(w(idx + 1:idx + H2), 1, H2);
    idx = idx + H2;

    % Extract b3 as scalar
    b3 = w(idx + 1);

    % Forward pass
    a1 = tanh(weights1 * inp + bias1);      % H1 x 1
    a2 = tanh(weights2 * a1 + b2);       % H2 x 1
    y  = W3 * a2 + b3;             % Scalar output
end



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
        %fprintf('No Detect \n');
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

function maxPd = simulateAndPlot(x, net)
    %% assemble the net
    % Assign weights/biases to a fresh copy of the NN
    nnFcn = @(n,Active) mySimpleNN(x, [n;Active], [3,3]);
    %net = setwb(net, x);

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
    Nd = 7.3327e+06;        % Number of detectors (must be odd)
    wd = 0.2750;               % Total detector array width [m]
    f = 0.0030;               % Focal length [m]
    delta_x_d = 3.7500e-08;       % Width of one detector element [m]
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

    % make missile run in random direction
    run_angle = rand()*2*pi;

    Rot = [cos(run_angle) -sin(run_angle); 
     sin(run_angle)  cos(run_angle)];

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
        %T_pos = [(-Vt_mag * t) + x0_t, y0_t];
        T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*sin(t/2) + .1*Vt_mag*sin(t*2) ) + y0_t];
        T_pos = [((-Vt_mag*t ) + x0_t )*cos(run_angle) - (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t*sin(run_angle), ((-Vt_mag*t ) + x0_t )*cos(run_angle) + (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t *sin(run_angle)];
        T_pos = (Rot * T_pos')';
        %T_pos_est = [(-Vt_mag * (t+1)) + x0_t, y0_t]; % consider using this to correct for laggining the target
        
        % Call Detect_IR to obtain the detector index and active pixel count.
        % Function call: [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, 
        %           r_tgt_missile, L_target, Vt_mag, x0_t, y0_t)
        [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, r_tgt_missile, L_target, Vt_mag, x0_t, y0_t);
    
        % --- 2. NN turn‑rate command
        inp = [n; Active];                   % 2×1
        %cmd_turn_rate = net(inp);                 % network output (1×1)
        cmd_turn_rate = nnFcn(n, Active);  % where nnFcn = @(n,Active) mySimpleNN(x, [n;Active], layerSizes)


        % --- 3. Limit and apply turn rate
        % Enforce the maximum turn rate limit.
        if abs(cmd_turn_rate) > max_turn_rate
            cmd_turn_rate = sign(cmd_turn_rate) * max_turn_rate;
        end
        
        % Update interceptor heading using the commanded turn rate.
        current_angle = atan2(V(2), V(1));       % Current interceptor heading [rad]
        new_angle = current_angle + cmd_turn_rate * dt;
        V = speed * [cos(new_angle), sin(new_angle)];
        
        % Update interceptor position using new velocity.
        S = S + V * dt;
        
        % Calculate fuze/damage probability using the Fuze function:
        % Function call: [Pd] = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag)
        Pd = Fuze(V, r_tgt_missile, L_target, K, S, phi, t, x0_t, y0_t, Vt_mag, T_pos);
        
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
            drawnow;
            % pause(0.1)
        end
    
        
        % Early exit on intercept (optional)
        if Pd > 0.9
            break;
        end
    
        % Check for collision: if the distance is less than the threshold, stop simulation.
        if norm(S - T_pos) < collision_threshold
            collisionDetected = true;
            disp('Contact detected! Exiting simulation.');
            break;  % Exit simulation loop.
        end
        
    
        % % Check for pause toggle
        while getappdata(hFig, 'isPaused')
            pause(0.1);
        end
    end

    maxPd = max(Pd_history);

    maxPd = -maxPd; % negative is good for minimization
    close
end