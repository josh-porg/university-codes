
% enablePlottting
enablePlottting = false; % should it p.lot as it goes?
add_noise = true; % add noise to sensor inputs for training?

%%%%%%%%%%
% TODO: used seeded random evasion. store the seed and evaluate the entire
% generation on the same seeds. regegerate seeds for next generation
%%%%%%%%%%

%%%%%%%%%%
%TODO: use len_hist = x(end); and make x one larger than needed
%%%%%%%%%%

% -------------------------------------------------------------------------
% 1. Create a feedforward network with two hidden layers of 32 neurons each
% -------------------------------------------------------------------------
inputsize = 8+20;
layerSizes = [2,2,2,2,2,2,2,2,2,2,2,2].^5;
outputsize = 20;

inputsize = 13;
layerSizes = [40,40,40,40,40];
layerSizes = [9,9,9];
%layerSizes = [15,15,15,15];

outputsize = 4;

if outputsize > 1
    inputsize = inputsize + 2*outputsize-1;
end

nVars = layerSizes(1)*inputsize + layerSizes(1);
for i = 2:length(layerSizes)
    layerSize = layerSizes(i);
    nVars = nVars + layerSizes(i)*layerSizes(i-1) + layerSizes(i); % weights + bias
end
nVars = nVars + outputsize*layerSizes(end) + outputsize; % weights + bias
nVars = nVars+1; % add one for the smoothing size and 7 of the other smoothing sizes

% nVars = layerSizes(1)*2 + layerSizes(1) ...     % weights1 + bias1
%       + layerSizes(2)*layerSizes(1) + layerSizes(2) ... % weights2 + b2
%       + layerSizes(2) + 1;                       % W3 + b3

% % 4. Set GA bounds (e.g. ±1 around zero)
lb = -1*ones(1, nVars);
ub =  1*ones(1, nVars);

% 5. Create a fitness-function handle (passing extra args)
fitnessFcn = @(x) computeFitness(x, enablePlottting, add_noise);

% 6. Configure GA options
popSize = 40; % how many individual per generation
elitePercent = .10; % save this fraction of population (the best)

opts = optimoptions('ga', ...
    'PopulationSize', popSize, ...
    'MaxGenerations', 70, ...
    'UseParallel', false ...
    , ...
    'EliteCount', ceil(elitePercent*popSize), ...
    'Display', 'iter');

% 7. Run GA!
[xBest, fvalBest] = ga(fitnessFcn, nVars, [], [], [], [], lb, ub, [], opts);

% 8. Build the optimized network

bestNN = @(n,Active) mySimpleNN(xBest, [n;Active], layerSizes);
% netOptimized = setwb(netTemplate, xBest);
% disp(['Optimized max Pd = ', num2str(1 - fvalBest)]);

maxPd = simulateAndPlot(xBest, true, false);

% plot a few to entertain
for i = 1:50
    maxPd = simulateAndPlot(xBest, true, false);
end


function Pd_avg = computeFitness(x, enablePlottting, add_noise)
    Pd_avg = 0;
    n_inters = 15;
    % do it 3 times because the running is random
    for itteration = 1:n_inters
        Pd_avg = Pd_avg + simulateAndPlot(x, enablePlottting, add_noise);
    end

    Pd_avg = Pd_avg/n_inters;
    close all;
end




%% mySimpleNN.m
% A feed‑forward network with two hidden layers you control by a weight vector.
% Usage: y = mySimpleNN(w, inp, layerSizes)
%
%   w          : 1×N vector of all weights & biases
%   inp        : 2×1 input vector [n; Active]
%   layerSizes : [H1, H2] sizes of hidden layers
%   y          : scalar network output

function y = mySimpleNN(w, inp, layerSizes, nOutputs)
    % Determine input size
    inputSize = length(inp);
    
    % Initialize index
    idx = 0;
    
    % Initialize layer outputs
    layerInput = inp;
    
    % Loop through each layer (except the last one)
    for i = 1:length(layerSizes) - 1
        % Current layer size (number of neurons)
        currentLayerSize = layerSizes(i);
        
        % Next layer size
        nextLayerSize = layerSizes(i + 1);
        
        % Extract weights and reshape
        nWeights = currentLayerSize * inputSize;
        weights = reshape(w(idx + 1 : idx + nWeights), currentLayerSize, inputSize);
        idx = idx + nWeights;
        
        % Extract bias and reshape
        bias = reshape(w(idx + 1 : idx + currentLayerSize), currentLayerSize, 1);
        idx = idx + currentLayerSize;
        
        % Forward pass (tanh activation)
        layerInput = tanh(weights * layerInput + bias);
        %layerInput = max(0,weights * layerInput + bias); % RELU activation
        %layerInput = sigmoid(weights * layerInput + bias);
        
        % Update input size for next layer
        inputSize = currentLayerSize;
    end
    
    % Final (output) layer (linear activation)
    % Now supports nOutputs instead of fixed 1 output
    W_final = reshape(w(idx + 1 : idx + layerSizes(end)*nOutputs), nOutputs, layerSizes(end));
    idx = idx + layerSizes(end)*nOutputs;
    b_final = reshape(w(idx + 1 : idx + nOutputs), nOutputs, 1);
    
    % Output (now nOutputs x 1 vector)
    y = W_final * layerInput + b_final;
end

function [y, lstmStates] = myLSTMNN(w, inpSequence, layerSizes, prevStates)
    % Inputs:
    %   w: Weight vector (flattened LSTM + FC weights)
    %   inpSequence: [inputSize x sequenceLength] matrix (time-series data)
    %   layerSizes: [inputSize, lstmHiddenUnits, fcLayer1, ..., outputSize]
    %   prevStates: Struct with fields {h, C} (hidden/cell states from last timestep)
    %
    % Outputs:
    %   y: Final output (same as mySimpleNN)
    %   lstmStates: Updated {h, C} for next call

    %% --- LSTM Layer ---
    inputSize = layerSizes(1);
    lstmHiddenUnits = layerSizes(2);
    
    % Extract LSTM weights from 'w'
    idx = 0;
    gateSize = lstmHiddenUnits * (inputSize + lstmHiddenUnits);
    
    W_f = reshape(w(idx + 1 : idx + gateSize), [lstmHiddenUnits, inputSize + lstmHiddenUnits]);
    idx = idx + gateSize;
    W_i = reshape(w(idx + 1 : idx + gateSize), [lstmHiddenUnits, inputSize + lstmHiddenUnits]);
    idx = idx + gateSize;
    W_C = reshape(w(idx + 1 : idx + gateSize), [lstmHiddenUnits, inputSize + lstmHiddenUnits]);
    idx = idx + gateSize;
    W_o = reshape(w(idx + 1 : idx + gateSize), [lstmHiddenUnits, inputSize + lstmHiddenUnits]);
    idx = idx + gateSize;
    
    % Biases
    b_f = reshape(w(idx + 1 : idx + lstmHiddenUnits), [lstmHiddenUnits, 1]);
    idx = idx + lstmHiddenUnits;
    b_i = reshape(w(idx + 1 : idx + lstmHiddenUnits), [lstmHiddenUnits, 1]);
    idx = idx + lstmHiddenUnits;
    b_C = reshape(w(idx + 1 : idx + lstmHiddenUnits), [lstmHiddenUnits, 1]);
    idx = idx + lstmHiddenUnits;
    b_o = reshape(w(idx + 1 : idx + lstmHiddenUnits), [lstmHiddenUnits, 1]);
    idx = idx + lstmHiddenUnits;
    
    % Initialize states if not provided
    if nargin < 4 || isempty(prevStates)
        h = zeros(lstmHiddenUnits, 1);
        C = zeros(lstmHiddenUnits, 1);
    else
        h = prevStates.h;
        C = prevStates.C;
    end
    
    % Process sequence through LSTM
    sequenceLength = size(inpSequence, 2);
    for t = 1:sequenceLength
        x_t = inpSequence(:, t);
        z = [h; x_t];  % Concatenate
        
        % LSTM gates
        f_t = sigmoid(W_f * z + b_f);
        i_t = sigmoid(W_i * z + b_i);
        C_tilde = tanh(W_C * z + b_C);
        o_t = sigmoid(W_o * z + b_o);
        
        % Update states
        C = f_t .* C + i_t .* C_tilde;
        h = o_t .* tanh(C);
    end
    
    % Save final states
    lstmStates.h = h;
    lstmStates.C = C;
    
    %% --- FC Layers (Same as mySimpleNN) ---
    layerInput = h;  % Output of LSTM -> input to FC layers
    inputSize = lstmHiddenUnits;
    
    for i = 3:length(layerSizes) - 1
        currentLayerSize = layerSizes(i);
        nextLayerSize = layerSizes(i + 1);
        
        % Extract weights and bias
        nWeights = currentLayerSize * inputSize;
        weights = reshape(w(idx + 1 : idx + nWeights), [currentLayerSize, inputSize]);
        idx = idx + nWeights;
        bias = reshape(w(idx + 1 : idx + currentLayerSize), [currentLayerSize, 1]);
        idx = idx + currentLayerSize;
        
        % Forward pass (tanh activation)
        layerInput = tanh(weights * layerInput + bias);
        inputSize = currentLayerSize;
    end
    
    % Output layer (linear)
    W_final = reshape(w(idx + 1 : idx + layerSizes(end)), [1, layerSizes(end)]);
    b_final = w(idx + layerSizes(end) + 1);
    y = W_final * layerInput + b_final;
end

% Helper: Sigmoid
function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end

function [n,Active] = Detect_IR(V,S,Nd,wd,f,delta_x_d,time_elapsed,...
    r_tgt_missile,L,Vt_mag,x0_t,y0_t, T_pos)
    % This function calculates the detector index (n) and the number of illuminated
    % detectors (Active) for the IR seeker using the current interceptor and target
    % parameters.
    
    % Compute target's current position relative to some fixed frame.
    St = T_pos;
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

function maxPd = simulateAndPlot(x, enablePlottting, add_noise)
    %% assemble the net
    % Assign weights/biases to a fresh copy of the NN
    nnFcn = @(inp) mySimpleNN(x, inp, [5,5,5,5,5,5]);
    nnFcn = @(inp) mySimpleNN(x, inp, [7,7,7,7,7,7]);
    nnFcn = @(inp) mySimpleNN(x, inp, [3,3,3]);
    layerShapes = [7,7,7,7,7,7]*2;
    layerShapes = [2,2,2,2,2,2,2,2,2,2,2,2].^5;
    layerShapes = [9,9,9];
    layerShapes = [15,15,15,15];
    layerShapes = [25,25,25,25];
    layerShapes = [25,25,25,25,25];
    layerShapes = [9,9,9];
    %layerShapes = [15,15,15,15];

    n_outputs = 4;
    output = zeros(n_outputs,1);

    %nnFN = @(inp) myLSTMNN(x, inp, [2, 64, 32, 1], prevStates)
    %net = setwb(net, x);

    %% Simulation Parameters
    dt = 0.1;              % Time step [s]
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
    x0_t = rand(1) * 2000;             % Target initial x-position [m]
    y0_t = rand(1) * 2000;               % Target y-position (constant) [m]
    Vt_mag = 200 + rand(1) *300;           % Target speed [m/s] (target velocity vector is [-Vt_mag, 0])
    
    %% IR Seeker Detector Parameters
    Nd = 7.3327e+06;        % Number of detectors (must be odd)
    wd = 0.2750;               % Total detector array width [m]
    f = 0.0030;               % Focal length [m]
    delta_x_d = 3.7500e-08;       % Width of one detector element [m]
    r_tgt_missile = 0.15;   % Radius of target missile [m] (30 cm diameter)
    L_target = 0.6;         % Length of target missile [m] (e.g., 60 cm)

    %% IR Seeker Detector Parameters
    Nd = 2.1998e+04;        % Number of detectors (must be odd)
    wd = 0.2750;               % Total detector array width [m]
    f = 0.0030;               % Focal length [m]
    delta_x_d = 1.2500e-05;       % Width of one detector element [m]
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

    % Initialize buffer (circular buffer for efficiency)
    buffer_size = 10;  % Store last 10 observations
    time_series_buffer = zeros(buffer_size, 2);  % Columns: [x_norm, y_norm, n_norm]
    
    %% Create Figure and Pause Button
    if enablePlottting
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
    end

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
    smoothNoise = smoothNoise / max(smoothNoise);
    
    %% stochastic evasion v2
    % Initialize the stochastic evation v2 parameters
    maxTurnRate = deg2rad(90); % Maximum turning rate (radians per second)
    maxAmplitude = 300; % Maximum deviation from base path

    % version 2:
    pathTypes = {"straight", "curved", "orbital", "sinusoidal", "squiggle", "random_turns"};
    weights = [0.15, 0.2, 0.05, 0.1, 0.2, 0.3]; % Probability distribution of path types
    weights = [0.15, 0.15, 0.05, 0.15, 0.15, 0.35]; % Probability distribution of path types
    weights = [0.2, 0.1, 0.05, 0.15, 0.15, 0.35]; % Probability distribution of path types
    weights = [0.25, 0.1, 0.05, 0.15, 0.15, 0.30]; % Probability distribution of path types

    pathType = pathTypes{find(rand < cumsum(weights), 1)};
    
    % Set parameters based on path type
    switch pathType
        case "straight"
            pathParams = struct('turnRate', 0, 'amplitude', 0);
        case "curved"
            pathParams = struct('turnRate', 0.5*maxTurnRate*(1-2*rand()), 'amplitude', 0);
        case "orbital"
            pathParams = struct('turnRate', 0.7*maxTurnRate + 0.3*maxTurnRate*(1-2*rand()), ...
                               'amplitude', NaN);
            pathParams.amplitude = 0.2*Vt_mag/pathParams.turnRate;
        case "sinusoidal"
            pathParams = struct('turnRate', 0.3*maxTurnRate*(1-2*rand()), ...
                              'amplitude', NaN);
            pathParams.amplitude = 0.3*Vt_mag/pathParams.turnRate;
        case "squiggle"
            pathParams = struct('turnRate', 0.3*maxTurnRate*(1-2*rand()), ...
                              'amplitude', NaN);
            pathParams.amplitude = 0.3*Vt_mag/pathParams.turnRate;
        case "random_turns"
            pathParams = struct('turnRate', 0.8*maxTurnRate, ... % Aggressive turns
                               'cruiseRate', 0.1*maxTurnRate*rand(), ... % Gentle cruise
                               'minCruiseTime', 1.0 + rand(), ... % Min time between turns [s]
                               'maxCruiseTime', 4.0 + 3*rand(), ... % Max time between turns [s]
                               'nextTurnTime', 0, ... % Will be initialized
                               'inTurnPhase', false, ... % Start in cruise
                               'turnDuration', 0.3 + 0.5*rand()); % How long turns last [s]
    % Set first turn time
    pathParams.nextTurnTime = pathParams.minCruiseTime + ...
                             (pathParams.maxCruiseTime-pathParams.minCruiseTime)*rand();
    end
    
    % Additional random parameters
    turnFrequency = 0.5 + rand(); % How often major turns occur
    baseAngle = run_angle; % Store initial angle

    

    % Initialize velocity and position
    currentVel = Vt_mag * [cos(run_angle); sin(run_angle)];
    currentPos = [x0_t; y0_t];

    % initialize commanded turn rate
    cmd_turn_rate = 0;

    % previous positions
    n_old = 0;
    Active_old = 0;

    %previous velocity
    dn_dt_old = 0;
    dActive_dt_old = 0;

    % integrals
    Active_int = 0;
    n_int = 0;
    LOS_int = 0;
    LOS_m_cmd_int = 0;

    % smoothing
    len_hist = 1 + floor(30 * (1 + x(end))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    LOS_hist = zeros(len_hist,1);
    LOS_rate_hist = zeros(len_hist,1);
    Active_hist = zeros(len_hist,1);
    Active_rate_hist = zeros(len_hist,1);
    cmd_turn_rate_hist = zeros(len_hist,1);
    current_angle_hist = zeros(len_hist,1);
    LOS_m_cmd_int_hist = zeros(len_hist,1);
    % len_hist1 = 1 + floor(30 * (1 + x(end-1))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    % len_hist2 = 1 + floor(30 * (1 + x(end-2))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    % len_hist3 = 1 + floor(30 * (1 + x(end-3))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    % len_hist4 = 1 + floor(30 * (1 + x(end-4))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    % len_hist5 = 1 + floor(30 * (1 + x(end-5))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    % len_hist6 = 1 + floor(30 * (1 + x(end-6))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    % len_hist7 = 1 + floor(30 * (1 + x(end-7))/2); % (1 + x(end))/2 is in range 0 to 1 thus this expresion is in tagenge 1 to 31
    % LOS_hist = zeros(len_hist1,1);
    % LOS_rate_hist = zeros(len_hist2,1);
    % Active_hist = zeros(len_hist3,1);
    % Active_rate_hist = zeros(len_hist4,1);
    % cmd_turn_rate_hist = zeros(len_hist5,1);
    % current_angle_hist = zeros(len_hist6,1);
    % LOS_m_cmd_int_hist = zeros(len_hist7,1);

    if n_outputs > 1
        output_hist = zeros(n_outputs, len_hist);
    end


    %% Simulation Loop
    for i = 1:numSteps
        t = time(i);

        % Update Target Position: Target moves leftward
        %T_pos = [(-Vt_mag * t) + x0_t, y0_t];
        T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*sin(t/2) + .1*Vt_mag*sin(t*2) ) + y0_t];
        T_pos = [((-Vt_mag*t ) + x0_t )*cos(run_angle) - (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t*sin(run_angle), ((-Vt_mag*t ) + x0_t )*cos(run_angle) + (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t *sin(run_angle)];
        T_pos = (Rot * T_pos')';
        %T_pos_est = [(-Vt_mag * (t+1)) + x0_t, y0_t]; % consider using this to correct for laggining the target

        %stochasic v2
        % Calculate desired velocity direction based on path type
        switch pathType
            case "straight"
                desiredDirection = run_angle;
            case "curved"
                % Gradual turn with limited rate
                currentVel = real(currentVel); % ensure real
                desiredDirection = atan2(currentVel(2), currentVel(1)) + pathParams.turnRate*dt;
            case "orbital"
                % Circular motion - direction always perpendicular to radius
                radius = currentPos - [x0_t; y0_t];
                if norm(radius) < eps
                    desiredDirection = run_angle;
                else
                    desiredDirection = atan2(radius(2), radius(1)) + pi/2;
                end
            case "sinusoidal"
                % Sinusoidal path along initial direction
                lateralOffset = pathParams.amplitude * sin(pathParams.turnRate * t);
                desiredDirection = run_angle + atan(lateralOffset/(Vt_mag*t));

                periodAdjustor1 = 1.5+rand(1);
                periodAdjustor2 = 1.5+rand(1);

                T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*sin(t/periodAdjustor1) + .1*Vt_mag*sin(t*periodAdjustor2) ) + y0_t];
                T_pos = (Rot * T_pos')';

            case "squiggle"
                % Sinusoidal path along initial direction
                lateralOffset = pathParams.amplitude * sin(pathParams.turnRate * t);
                desiredDirection = run_angle + atan(lateralOffset/(Vt_mag*t));

                T_pos = [((-Vt_mag*t ) + x0_t )*cos(run_angle) - (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t*sin(run_angle), ((-Vt_mag*t ) + x0_t )*cos(run_angle) + (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t *sin(run_angle)];
                T_pos = (Rot * T_pos')';
            case "random_turns"
                % Check if we should switch phases
                if pathParams.inTurnPhase
                    % Check if turn duration has elapsed
                    if (t - pathParams.turnStartTime) >= pathParams.turnDuration
                        pathParams.inTurnPhase = false;
                        pathParams.nextTurnTime = t + pathParams.minCruiseTime + ...
                                                (pathParams.maxCruiseTime-pathParams.minCruiseTime)*rand();
                    end
                else
                    % Check if it's time for a new turn
                    if t >= pathParams.nextTurnTime
                        pathParams.inTurnPhase = true;
                        pathParams.turnStartTime = t;
                        % Randomly choose turn direction (left or right)
                        pathParams.currentTurnSign = sign(randn());
                    end
                end
                
                % Calculate desired direction based on current phase
                if pathParams.inTurnPhase
                    % In turn phase - apply aggressive turn
                    currentVel = real(currentVel); % ensure real
                    desiredDirection = atan2(currentVel(2), currentVel(1)) + ...
                                      pathParams.currentTurnSign * pathParams.turnRate * dt;
                else
                    currentVel = real(currentVel); % ensure real
                    % In cruise phase - gentle curve or straight
                    desiredDirection = atan2(currentVel(2), currentVel(1)) + ...
                                      pathParams.cruiseRate * dt;
                end
        end
        
        % Current direction
        currentVel = real(currentVel); % ensure real
        currentDirection = atan2(currentVel(2), currentVel(1));
        
        % Calculate maximum allowed direction change
        maxDirectionChange = maxTurnRate * dt;
        desiredDirectionChange = desiredDirection - currentDirection;
        
        % Normalize angle difference to [-pi, pi]
        desiredDirectionChange = mod(desiredDirectionChange + pi, 2*pi) - pi;
        
        % Apply turn rate constraint
        actualDirectionChange = sign(desiredDirectionChange) * min(abs(desiredDirectionChange), maxDirectionChange);
        
        % Calculate new velocity direction
        newDirection = currentDirection + actualDirectionChange;
        
        % Update velocity vector (maintaining constant magnitude)
        currentVel = Vt_mag * [cos(newDirection); sin(newDirection)];
        currentVel = real(currentVel); % ensure real
        
        % Add evasive wiggles (perpendicular to current direction)
        wiggleMagnitude = 0.1 * Vt_mag;
        wiggleDirection = newDirection + pi/2; % Perpendicular wiggle
        wiggle = wiggleMagnitude * smoothNoise(i) * [cos(wiggleDirection); sin(wiggleDirection)];
        
        % Update position with constrained velocity + wiggles
        currentPos = currentPos + (currentVel + wiggle) * dt;
        
        % Final position output
        T_pos = currentPos';
        
        % Apply additional rotation if needed (your original Rot matrix)
        T_pos = (Rot * T_pos')';
        %end stocahstic v2

        if pathType == "sinusoidal"
            T_pos = [(-Vt_mag*t ) + x0_t, (Vt_mag*sin(t/2) + .1*Vt_mag*sin(t*2) ) + y0_t];
            T_pos = (Rot * T_pos')';
        end
        if pathType == "squiggle"
            T_pos = [((-Vt_mag*t ) + x0_t )*cos(run_angle) - (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t*sin(run_angle), ((-Vt_mag*t ) + x0_t )*cos(run_angle) + (Vt_mag*smoothNoise(i) + .1*Vt_mag*sin(t*2) ) + y0_t *sin(run_angle)];
            T_pos = (Rot * T_pos')';
        end
        
        % Call Detect_IR to obtain the detector index and active pixel count.
        % Function call: [n, Active] = Detect_IR(V, S, Nd, wd, f, delta_x_d, t, 
        %           r_tgt_missile, L_target, Vt_mag, x0_t, y0_t)
        %[n, Active] = Detect_IR_v3(V, S, Nd, wd, f, delta_x_d, t, r_tgt_missile, L_target, Vt_mag, x0_t, y0_t, T_pos);
        [n, Active] = Detect_IR_v3(V, S, Nd, wd, f, delta_x_d, t, r_tgt_missile, L_target, Vt_mag, x0_t, y0_t, T_pos);

        % add Noise (might improve training)
        if add_noise
            n = n + randn() * .01 * Nd;
            Active = max(Active + randn() * .01 * Nd, 0);
        end
    
        % compute derivatives with backward difference
        dn_dt = (n - n_old) / dt;
        dActive_dt = (Active - Active_old) / dt;

        d2ndt2 = (dn_dt - dn_dt_old) / dt;
        d2Active_dt2 = (dActive_dt - dActive_dt_old) / dt;

        % integrals
        Active_int = Active_int + Active/Nd * dt;
        n_int = n_int + n * dt;

        LOS = -atan2(n * delta_x_d, f);
        LOS_rate = (-atan2(n * delta_x_d, f) - -atan2(n_old * delta_x_d, f)) / dt;
        LOS_int = LOS_int + -atan2(n * delta_x_d, f) * dt;

        LOS_m_cmd_int = LOS_m_cmd_int + (LOS - cmd_turn_rate)/dt;

        V = real(V); % keep it real
        current_angle = atan2(V(2), V(1));       % Current interceptor heading [rad]

        % update the histories
        LOS_hist = [LOS; LOS_hist(1:end-1)];
        LOS_rate_hist = [LOS_rate; LOS_rate_hist(1:end-1)];
        LOS_m_cmd_int_hist = [LOS_m_cmd_int; LOS_m_cmd_int_hist(1:end-1)];
        cmd_turn_rate_hist = [cmd_turn_rate; cmd_turn_rate_hist(1:end-1)];
        Active_hist = [Active/Nd; Active_hist(1:end-1)];
        Active_rate_hist = [dActive_dt/Nd; Active_rate_hist(1:end-1)];
        current_angle_hist = [current_angle/(2*pi); current_angle_hist(1:end-1)];

        if n_outputs > 1
            output_hist = [output, output_hist(:, 1:end-1)];
        end

        % smoothed parameters
        LOS_smooth = mean(LOS_hist);
        LOS_rate_smooth = mean(LOS_rate_hist);
        LOS_m_cmd_int_smooth = mean(LOS_m_cmd_int_hist);
        cmd_turn_rate_smooth = mean(cmd_turn_rate_hist);
        Active_smooth = mean(Active_hist);
        Active_rate_smooth = mean(Active_rate_hist);
        current_angle_smooth = mean(current_angle_hist);
        
        if n_outputs > 1
            %smooth output
            output_smooth = mean(output_hist,2);

            % differentialte output then smooth it
            output_rate_smooth = mean(diff(output_hist)/dt,2);
        end


        

        % % Update buffer (timestep history of the missile)
        % current_observation = [n/Nd, Active/Nd];
        % time_series_buffer = [current_observation; time_series_buffer(1:end-1, :)];  % Push new, pop oldest
        % 
        % % Calculate velocity (Δx, Δy) between consecutive steps
        % delta_x = real(diff(time_series_buffer(:, 1))/dt);  
        % delta_y = real(diff(time_series_buffer(:, 2))/dt);  
        % 
        % % Pad with zeros to maintain buffer size
        % delta_x = [0; delta_x];  % No delta for the first observation  
        % delta_y = [0; delta_y];  
        % 
        % % Augment buffer with deltas
        % augmented_buffer = [time_series_buffer, delta_x, delta_y];  
        
        % --- 2. NN turn‑rate command
        % inp = [n; Active];                   % 2×1 - works well with a 3x3 network
        % inp = [-atan2(n * delta_x_d, f); Active/Nd]; % this also works well
        % 
        % inp = [-atan2(n * delta_x_d, f); Active/Nd; (-atan2(n * delta_x_d, f) - -atan2(n_old * delta_x_d, f)) / dt; (Active/Nd - Active_old/Nd) / dt]; % this is one of the bestter inputs (run using [7,7,7,7,7,7] or [7,7,7,7,7,7]*2 with 7 outputs )
        % inp = [-atan2(n * delta_x_d, f); Active/Nd; (-atan2(n * delta_x_d, f) - -atan2(n_old * delta_x_d, f)) / dt; (Active/Nd - Active_old/Nd) / dt; current_angle/(2*pi); cmd_turn_rate/deg2rad(max_turn_rate)]; % this is one of the bestter inputs (run using [7,7,7,7,7,7]*2 with 7 outputs)
        % inp = [-atan2(n * delta_x_d, f); Active/Nd; (-atan2(n * delta_x_d, f) - -atan2(n_old * delta_x_d, f)) / dt; (Active/Nd - Active_old/Nd) / dt; current_angle/(2*pi); cmd_turn_rate/deg2rad(max_turn_rate); t/60]; % this is one of the bestter inputs
        % inp = [-atan2(n * delta_x_d, f); Active/Nd; (-atan2(n * delta_x_d, f) - -atan2(n_old * delta_x_d, f)) / dt; (Active/Nd - Active_old/Nd) / dt; current_angle/(2*pi); cmd_turn_rate/deg2rad(max_turn_rate); LOS_int]; % this is one of the bestter inputs (run using [7,7,7,7,7,7]*2 with 7 outputs)
        % 
        % inp = [LOS; LOS_smooth; LOS_rate; LOS_rate_smooth; LOS_int; Active/Nd; Active_smooth; (Active/Nd - Active_old/Nd) / dt; Active_rate_smooth; current_angle/(2*pi); current_angle_smooth; cmd_turn_rate/deg2rad(max_turn_rate); cmd_turn_rate_smooth; output_rate_smooth];
        % inp = [LOS; LOS_smooth; LOS_rate; LOS_rate_smooth; LOS_int; Active/Nd; Active_smooth; (Active/Nd - Active_old/Nd) / dt; Active_rate_smooth; current_angle/(2*pi); current_angle_smooth; cmd_turn_rate/deg2rad(max_turn_rate); cmd_turn_rate_smooth; t/60];
        % inp = [LOS; LOS_smooth; LOS_rate; LOS_rate_smooth; LOS_int; Active/Nd; Active_smooth; (Active/Nd - Active_old/Nd) / dt; Active_rate_smooth; current_angle/(2*pi); current_angle_smooth; cmd_turn_rate/deg2rad(max_turn_rate); cmd_turn_rate_smooth;];
        % inp = [LOS; LOS_smooth; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; LOS_int; Active/Nd; Active_smooth; Active_rate_smooth];

        %inp = [LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; Active/Nd; Active_smooth; Active_rate_smooth]; %good for porportional with layerShapes = [9,9,9]; 4 outupt returns
        %inp = [LOS; LOS_smooth; LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; Active/Nd; Active_smooth; Active_rate_smooth];
        % inp = [LOS; LOS_smooth; LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; LOS_m_cmd_int; LOS_m_cmd_int_smooth; Active/Nd; Active_smooth; Active_rate_smooth];
        %inp = [LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; LOS_m_cmd_int; LOS_m_cmd_int_smooth; Active/Nd; Active_smooth; Active_rate_smooth];
        inp = [(LOS - current_angle)/(2*pi); LOS_smooth/(2*pi) - current_angle_smooth; LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; LOS_m_cmd_int; LOS_m_cmd_int_smooth; Active/Nd; Active_smooth; Active_rate_smooth];

        % attempted proportional PID
        inp = [(LOS - current_angle)/(2*pi); LOS_smooth/(2*pi) - current_angle_smooth; LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; LOS_m_cmd_int; LOS_m_cmd_int_smooth; Active/Nd; Active_smooth; Active_rate_smooth];%currently kind of workds with layerSizes = [15,15,15,15]; outputsize = 4;

        % good for porportional with layerShapes = [9,9,9]; with 4 outputs
        % inp = [LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; Active/Nd; Active_smooth; Active_rate_smooth]; %good for porportional with layerShapes = [9,9,9]; 4 outupt returns
        %inp = [LOS; LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; Active/Nd; Active_smooth; Active_rate_smooth];
        inp = [LOS_rate; LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; Active/Nd; Active_smooth; Active_rate_smooth]; %good for porportional with layerShapes = [9,9,9]; 4 outupt returns

        %inp = [n; Active; dn_dt; dActive_dt]; % does work well (with a 5,5,5,5,5,5 net)
        %inp = [n; Active; dn_dt; dActive_dt; V'; cmd_turn_rate];
        if output > 1
            %inp = [inp; output(2:end); output_smooth; output_rate_smooth]; % does not work well (with a 5,5,5,5,5,5 net)
            inp = [inp; output(2:end); output_smooth]; % does not work well (with a 5,5,5,5,5,5 net)
        end
        %inp = [n; Active; dn_dt; dActive_dt; d2ndt2; d2Active_dt2]; % does not work well (with a 7,7,7,7,7,7 net)

        

        %cmd_turn_rate = net(inp);                 % network output (1×1)
        %cmd_turn_rate = nnFcn(inp);  % where nnFcn = @(inp) mySimpleNN(x, inp, layerSizes)
        output = mySimpleNN(x, inp, layerShapes, n_outputs);
        cmd_turn_rate = output(1);


        % --- 3. Limit and apply turn rate
        % Enforce the maximum turn rate limit.
        if abs(cmd_turn_rate) > max_turn_rate
            cmd_turn_rate = sign(cmd_turn_rate) * max_turn_rate;
        end
        
        % Update interceptor heading using the commanded turn rate.
        V = real(V); % keep it real
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

        % store the current values of n and Active as the old values (update
        % them 
        n_old = n;
        Active_old = Active;
        
        if (enablePlottting) && (mod(i,50) == 0)
            % Update plot for interceptor and target paths and current positions
            set(hPlot, 'XData', S_history(1:i,1), 'YData', S_history(1:i,2));
            set(hTarget, 'XData', T_history(1:i,1), 'YData', T_history(1:i,2));
            set(hInterceptor, 'XData', S(1), 'YData', S(2));
            set(hTargetPos, 'XData', T_pos(1), 'YData', T_pos(2));
            drawnow;
            % pause(0.1)
        end
    
        
        % Early exit on intercept (optional)
        if Pd > 0.999
            break;
        end
    
        % Check for collision: if the distance is less than the threshold, stop simulation.
        if norm(S - T_pos) < collision_threshold
            collisionDetected = true;
            disp('Contact detected! Exiting simulation.');
            break;  % Exit simulation loop.
        end
        
        if enablePlottting
            % % Check for pause toggle
            while getappdata(hFig, 'isPaused')
                pause(0.1);
            end
        end
    end

    maxPd = max(Pd_history);

    maxPd = -maxPd; % negative is good for minimization
    close
end