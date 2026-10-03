function Pd = simAndPlot_requested(x, enablePlottting, add_noise, save_data)
    %% assemble the net
    % Assign weights/biases to a fresh copy of the NN
    layerShapes = [9,9,9];
    %layerShapes = [15,15,15,15];

    n_outputs = 4;
    output = zeros(n_outputs,1);

    %nnFN = @(inp) myLSTMNN(x, inp, [2, 64, 32, 1], prevStates)
    %net = setwb(net, x);

    %% Simulation Parameters
    overspd_ratio = 10; % number of timesteps between guidance updates
    guidance_update_rate = 0.1; % Time step [s]
    dt = guidance_update_rate / overspd_ratio; % Time step [s]
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
    x0_t = 300;
    y0_t = 0;
    Vt_mag = 200;           % Target speed [m/s] (target velocity vector is [-Vt_mag, 0])

    %% IR Seeker Detector Parameters
    Nd = 2.1998e+04;        % Number of detectors (must be odd)
    wd = 0.2750;               % Total detector array width [m]
    f = 0.0030;               % Focal length [m]
    delta_x_d = 1.2500e-05;       % Width of one detector element [m]
    r_tgt_missile = 0.15;   % Radius of target missile [m] (30 cm diameter)
    L_target = 0.6;         % Length of target missile [m] (e.g., 60 cm)
    
    %% more guidance and fuzing parameters
    collision_threshold = r_tgt_missile;
    
    %% Fuze (Damage) Parameters
    K = 90;                 % Number of fragments
    %phi = deg2rad(30);      % Fragmentation half-angle [rad] (example value)

    phi = deg2rad(15 + floor(45 * (1 + x(end-1))/2)); % let it choose its own phi ( Fragmentation half-angle [rad] )
    
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

    weights = [1, 0, 0, 0, 0, 0]; % Probability distribution of path types

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
            % n = n + randn() * .01 * Nd;
            % Active = max(Active + randn() * .01 * Nd, 0);
            n = n + randn() * .005 * Nd;
            Active = max(Active + randn() * .005 * Nd, 0);
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

        if mod(i,overspd_ratio) == 0
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
            %inp = [LOS_rate-cmd_turn_rate; LOS_rate_smooth-cmd_turn_rate; Active/Nd; Active_smooth; Active_rate_smooth]; %good for porportional with layerShapes = [9,9,9]; 4 outupt returns

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
        end


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
    
        if estimatePd(Active, n, rad2deg(phi), delta_x_d, f, r_tgt_missile, L_target, K) >= 0.9
            break;
        end
        
        if enablePlottting
            % % Check for pause toggle
            while getappdata(hFig, 'isPaused')
                pause(0.1);
            end
        end
    end

    Pd = Pd_history(i);
    close

    %% Data Saving Logic
    if save_data
        filename = 'simulation_data.mat';
        
        % Check if file exists and load existing data
        if exist(filename, 'file')
            loaded_data = load(filename);
            all_n_history = loaded_data.all_n_history;
            all_Active_history = loaded_data.all_Active_history;
            all_Pd_history = loaded_data.all_Pd_history;
            all_Phi_history = loaded_data.all_Phi_history;
        else
            all_n_history = [];
            all_Active_history = [];
            all_Pd_history = [];
            all_Phi_history = [];
        end
        
        % Truncate to actual simulation steps (in case of early exit)
        actual_steps = i;
        
        % Append new data
        all_n_history = [all_n_history; n_history(1:actual_steps)];
        all_Active_history = [all_Active_history; Active_history(1:actual_steps)];
        all_Pd_history = [all_Pd_history; Pd_history(1:actual_steps)];
        all_Phi_history = [all_Phi_history; phi* ones(size(Pd_history(1:actual_steps)))];
        
        % Save combined data
        save(filename, 'all_n_history', 'all_Active_history', 'all_Pd_history', 'all_Phi_history');
        
        % Display stats
        %fprintf('Saved data. Total samples: %d\n', length(all_n_history));
    end
end