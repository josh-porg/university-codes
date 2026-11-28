%% Custom RL Environment
classdef MissileGuidanceEnv < rl.env.MATLABEnvironment
    properties
        ResetFcn="reset"

        % Simulation parameters (same as original)
        enablePlotting = false;
        add_noise = true;
        overspd_ratio = 1;%10;
        guidance_update_rate = 0.1;
        dt = 0.1 / 10;
        T_total = 60;
        speed = 400;
        max_turn_deg = 100;
        max_turn_rate = deg2rad(100);
        r_tgt_missile = 0.15;
        L_target = 0.6;
        K = 90;
        collision_threshold = 0.15;
        
        % Neural network parameters
        inputSize = 10;
        layerSizes = [25, 25, 25];
        outputSize = 4;
        
        % State variables
        S = [0, 0];
        V = [0, 400];
        T_pos = [0, 0];
        x0_t = 0;
        y0_t = 0;
        Vt_mag = 200 + rand(); %* 300; % initializing this to 0 breaks things
        run_angle = 0;
        n_old = 0;
        Active_old = 0;
        cmd_turn_rate = 0;
        time = 0;
        stepCount = 0;
        phi = deg2rad(30);

        %integral
        LOS_m_cmd_int = 0;
        
        % History buffers
        len_hist = 10;
        LOS_hist = zeros(10,1);
        LOS_rate_hist = zeros(10,1);
        Active_hist = zeros(10,1);
        Active_rate_hist = zeros(10,1);
        cmd_turn_rate_hist = zeros(10,1);
        current_angle_hist = zeros(10,1);
        LOS_m_cmd_int_hist = zeros(10,1);
        
        % Target evasion parameters
        pathType = "straight";
        smoothNoise = [];
        pathParams = struct('turnRate', 0, 'amplitude', 0);
    end
    
    methods
        function this = MissileGuidanceEnv(enablePlotting, add_noise, inputSize, layerSizes, outputSize)

            % Define observation and action info
            obsInfo = rlNumericSpec([inputSize 1]);
            actInfo = rlNumericSpec([outputSize 1], 'LowerLimit', -1, 'UpperLimit', 1);

            % Initialize environment with observation and action info
            this = this@rl.env.MATLABEnvironment(obsInfo, actInfo);
            this.ResetFcn = @reset;

            % Initialize environment
            %this = this@rl.env.MATLABEnvironment();
            
            % Set properties
            this.enablePlotting = enablePlotting;
            this.add_noise = add_noise;
            this.inputSize = inputSize;
            this.layerSizes = layerSizes;
            this.outputSize = outputSize;
            
            % Define observation and action info
            obsInfo = rlNumericSpec([inputSize 1]);
            actInfo = rlNumericSpec([outputSize 1], 'LowerLimit', -1, 'UpperLimit', 1);
            
            this.ActionInfo = actInfo;
            this.ObservationInfo = obsInfo;
            
            % Initialize target evasion noise
            this.initializeNoise();
        end
        
        function initializeNoise(this)
            % Generate smooth noise for target evasion
            n = length(0:this.dt:this.T_total);
            cutoffFreq = .001;
            noise = randn(1, n);
            [b, a] = butter(2, cutoffFreq, 'low');
            this.smoothNoise = filtfilt(b, a, noise);
            this.smoothNoise = this.smoothNoise / max(this.smoothNoise);
        end
        
        function [obs, reward, done, info] = step(this, action)
            % Execute one time step in the environment
            
            % Update time
            this.time = this.time + this.dt;
            this.stepCount = this.stepCount + 1;
            
            % Update target position (stochastic evasion)
            this.updateTargetPosition();
            
            % Simulate IR detection
            [n, Active] = this.detectIR();
            
            % Add noise if enabled
            if this.add_noise
                n = n + randn() * 0.005 * this.getNd();
                Active = max(Active + randn() * 0.005 * this.getNd(), 0);
            end
            
            % Compute derivatives
            dn_dt = (n - this.n_old) / this.dt;
            dActive_dt = (Active - this.Active_old) / this.dt;
            
            % Compute LOS and rates
            LOS = -atan2(n * this.getDeltaXd(), this.getFocalLength());
            LOS_rate = (-atan2(n * this.getDeltaXd(), this.getFocalLength()) - ...
                       -atan2(this.n_old * this.getDeltaXd(), this.getFocalLength())) / this.dt;
            
            % Update integrals
            this.LOS_m_cmd_int = this.LOS_m_cmd_int + (LOS - this.cmd_turn_rate)/this.dt;
            
            % Get current angle
            current_angle = atan2(this.V(2), this.V(1));
            
            % Update history buffers
            if mod(this.stepCount, this.overspd_ratio) == 0
                this.updateHistoryBuffers(LOS, LOS_rate, Active, current_angle);
                
                % Get smoothed values
                LOS_smooth = mean(this.LOS_hist);
                LOS_rate_smooth = mean(this.LOS_rate_hist);
                LOS_m_cmd_int_smooth = mean(this.LOS_m_cmd_int_hist);
                cmd_turn_rate_smooth = mean(this.cmd_turn_rate_hist);
                Active_smooth = mean(this.Active_hist);
                Active_rate_smooth = mean(this.Active_rate_hist);
                current_angle_smooth = mean(this.current_angle_hist);
                
                % Create observation vector
                obs = [
                    %(LOS - current_angle)/(2*pi);
                    %LOS_smooth/(2*pi) - current_angle_smooth;
                    LOS_rate;
                    LOS_rate - this.cmd_turn_rate;
                    LOS_rate_smooth - this.cmd_turn_rate;
                    %this.LOS_m_cmd_int;
                    %LOS_m_cmd_int_smooth;
                    Active/this.getNd();
                    Active_smooth;
                    Active_rate_smooth
                ];
            else
                obs = this.getLastObservation();
            end
            

            if isa(action(1), "cell")
                action = cell2mat(action);
            end

            % Apply action (turn rate command)
            this.cmd_turn_rate = action(1) * this.max_turn_rate;
            
            % Limit turn rate
            if abs(this.cmd_turn_rate) > this.max_turn_rate
                this.cmd_turn_rate = sign(this.cmd_turn_rate) * this.max_turn_rate;
            end
            
            % Update interceptor heading and position
            new_angle = current_angle + this.cmd_turn_rate * this.dt;
            this.V = this.speed * [cos(new_angle), sin(new_angle)];
            this.S = this.S + this.V * this.dt;
            
            % Calculate damage probability
            Pd = this.calculateDamageProbability();
            
            % Store current values
            this.n_old = n;
            this.Active_old = Active;
            
            % Check termination conditions
            done = this.time >= this.T_total || Pd > 0.99 || ...
                   norm(this.S - this.T_pos) < this.collision_threshold;

            if done
                fprintf("Done: time: %g; ", this.time)
                fprintf("Pd: %g; ", Pd)
                fprintf("dist: %g\n", norm(this.S - this.T_pos))
            end

            
            % Reward is negative of max Pd (we want to maximize Pd)
            reward = -Pd;

            % (2) Improved reward shaping
            distance = norm(this.S - this.T_pos);
            reward = -Pd - 0.1 / (distance + 1e-5); % Bonus for proximity
            reward = Pd;  % Now reward ∈ [0, 1] (directly maximize Pd)
            % OR
            %reward = 10 * Pd - 1;  % Reward ∈ [-1, 9] (centered around 0)
            % New reward: Directly maximize Pd + distance bonus
            reward = Pd + 0.1 / (distance + 1e-5);  % ∈ [0, 1 + bonus]
            
            % Additional info
            info = struct('Pd', Pd, 'time', this.time);

            if done
                this.reset(); % manually rest
            end
        end
        
        function [obs, info] = reset(this)
            % Reset environment to initial state
            
            % Reset time and step count
            this.time = 0;
            this.stepCount = 0;
            
            % Reset interceptor state
            init_heading = pi/2;
            this.S = [0, 0];
            this.V = this.speed * [cos(init_heading), sin(init_heading)];
            
            % Reset target state with random parameters
            this.run_angle = rand() * 2 * pi;
            this.x0_t = rand() * 2000;
            this.y0_t = rand() * 2000;
            this.Vt_mag = 200 + rand() * 300;
            
            % Reset target evasion path
            this.initializeTargetPath();
            
            % Reset detection variables
            this.n_old = 0;
            this.Active_old = 0;
            this.cmd_turn_rate = 0;

            % reset the integral
            this.LOS_m_cmd_int = 0;
            
            % Reset history buffers
            this.LOS_hist = zeros(this.len_hist,1);
            this.LOS_rate_hist = zeros(this.len_hist,1);
            this.Active_hist = zeros(this.len_hist,1);
            this.Active_rate_hist = zeros(this.len_hist,1);
            this.cmd_turn_rate_hist = zeros(this.len_hist,1);
            this.current_angle_hist = zeros(this.len_hist,1);
            this.LOS_m_cmd_int_hist = zeros(this.len_hist,1);

            % fprintf("New episode: time: %g; ", this.time)
            %     %fprintf("Pd: %g; ", Pd)
            % fprintf("dist: %g; ", norm(this.S - this.T_pos))
            % fprintf('Target=(%.1f, %.1f), Speed=%.1f m/s, Angle=%.1f rad\n', ...
            % this.x0_t, this.y0_t, this.Vt_mag, this.run_angle);
            
            % Get initial observation
            obs = this.getInitialObservation();
            info = NaN;
        end
        
        function initializeTargetPath(this)
            % Initialize target evasion path type and parameters
            pathTypes = ["straight", "curved", "orbital", "sinusoidal", "squiggle", "random_turns"];
            weights = [0.25, 0.1, 0.05, 0.15, 0.15, 0.30];
            this.pathType = pathTypes{find(rand() < cumsum(weights), 1)};
            
            maxTurnRate = deg2rad(90);
            
            switch this.pathType
                case "straight"
                    this.pathParams = struct('turnRate', 0, 'amplitude', 0);
                case "curved"
                    this.pathParams = struct('turnRate', 0.5*maxTurnRate*(1-2*rand()), 'amplitude', 0);
                case "orbital"
                    this.pathParams = struct('turnRate', 0.7*maxTurnRate + 0.3*maxTurnRate*(1-2*rand()), ...
                                           'amplitude', 0);
                    this.pathParams = struct('turnRate', 0.7*maxTurnRate + 0.3*maxTurnRate*(1-2*rand()), ...
                                           'amplitude', 0.2*this.Vt_mag/this.pathParams.turnRate);
                case "sinusoidal"
                    this.pathParams = struct('turnRate', 0.3*maxTurnRate*(1-2*rand()), ...
                                           'amplitude', 0);
                    this.pathParams = struct('turnRate', 0.3*maxTurnRate*(1-2*rand()), ...
                                           'amplitude', 0.3*this.Vt_mag/this.pathParams.turnRate);
                case "squiggle"
                    this.pathParams = struct('turnRate', 0.3*maxTurnRate*(1-2*rand()), ...
                                           'amplitude', 0);
                    this.pathParams = struct('turnRate', 0.3*maxTurnRate*(1-2*rand()), ...
                                           'amplitude', 0.3*this.Vt_mag/this.pathParams.turnRate);
                case "random_turns"
                    this.pathParams = struct('turnRate', 0.3*maxTurnRate*(1-2*rand()), ...
                                           'amplitude', 0);
                    this.pathParams = struct('turnRate', 0.8*maxTurnRate, ...
                                           'cruiseRate', 0.1*maxTurnRate*rand(), ...
                                           'minCruiseTime', 1.0 + rand(), ...
                                           'maxCruiseTime', 4.0 + 3*rand(), ...
                                           'nextTurnTime', 1.0 + rand()*3, ...
                                           'inTurnPhase', false, ...
                                           'turnDuration', 0.3 + 0.5*rand());
            end
        end
        
        function updateTargetPosition(this)
            % Update target position based on evasion strategy
            
            % Calculate desired velocity direction based on path type
            switch this.pathType
                case "straight"
                    desiredDirection = this.run_angle;
                case "curved"
                    currentVel = this.Vt_mag * [cos(this.run_angle); sin(this.run_angle)];
                    desiredDirection = atan2(currentVel(2), currentVel(1)) + this.pathParams.turnRate*this.dt;
                case "orbital"
                    radius = [this.T_pos(1) - this.x0_t; this.T_pos(2) - this.y0_t];
                    if norm(radius) < eps
                        desiredDirection = this.run_angle;
                    else
                        desiredDirection = atan2(radius(2), radius(1)) + pi/2;
                    end
                case "sinusoidal"
                    periodAdjustor1 = 1.5 + rand();
                    periodAdjustor2 = 1.5 + rand();
                    this.T_pos = [(-this.Vt_mag*this.time) + this.x0_t, ...
                                 (this.Vt_mag*sin(this.time/periodAdjustor1) + 0.1*this.Vt_mag*sin(this.time*periodAdjustor2)) + this.y0_t];
                    Rot = [cos(this.run_angle) -sin(this.run_angle); 
                           sin(this.run_angle)  cos(this.run_angle)];
                    this.T_pos = (Rot * this.T_pos')';
                    return;
                case "squiggle"
                    this.T_pos = [((-this.Vt_mag*this.time) + this.x0_t)*cos(this.run_angle) - ...
                                 (this.Vt_mag*this.smoothNoise(this.stepCount) + 0.1*this.Vt_mag*sin(this.time*2)) + this.y0_t*sin(this.run_angle), ...
                                 ((-this.Vt_mag*this.time) + this.x0_t)*cos(this.run_angle) + ...
                                 (this.Vt_mag*this.smoothNoise(this.stepCount) + 0.1*this.Vt_mag*sin(this.time*2)) + this.y0_t*sin(this.run_angle)];
                    Rot = [cos(this.run_angle) -sin(this.run_angle); 
                           sin(this.run_angle)  cos(this.run_angle)];
                    this.T_pos = (Rot * this.T_pos')';
                    return;
                case "random_turns"
                    % Check if we should switch phases
                    if this.pathParams.inTurnPhase
                        % Check if turn duration has elapsed
                        if (this.time - this.pathParams.turnStartTime) >= this.pathParams.turnDuration
                            this.pathParams.inTurnPhase = false;
                            this.pathParams.nextTurnTime = this.time + this.pathParams.minCruiseTime + ...
                                                         (this.pathParams.maxCruiseTime-this.pathParams.minCruiseTime)*rand();
                        end
                    else
                        % Check if it's time for a new turn
                        if this.time >= this.pathParams.nextTurnTime
                            this.pathParams.inTurnPhase = true;
                            this.pathParams.turnStartTime = this.time;
                            this.pathParams.currentTurnSign = sign(randn());
                        end
                    end
                    
                    % Calculate desired direction based on current phase
                    if this.pathParams.inTurnPhase
                        % In turn phase - apply aggressive turn
                        desiredDirection = atan2(this.V(2), this.V(1)) + ...
                                          this.pathParams.currentTurnSign * this.pathParams.turnRate * this.dt;
                    else
                        % In cruise phase - gentle curve or straight
                        desiredDirection = atan2(this.V(2), this.V(1)) + ...
                                          this.pathParams.cruiseRate * this.dt;
                    end
            end
            
            % Update target position
            this.T_pos = this.T_pos + this.Vt_mag * [cos(desiredDirection), sin(desiredDirection)] * this.dt;
        end
        
        function [n, Active] = detectIR(this)
            % Simulate IR detection (simplified version)
            % Calculate relative position
            relPos = this.T_pos - this.S;
            
            % Calculate angle to target
            angleToTarget = atan2(relPos(2), relPos(1)) - atan2(this.V(2), this.V(1));
            
            % Simulate detector array
            Nd = this.getNd();
            delta_x_d = this.getDeltaXd();
            f = this.getFocalLength();
            
            % Calculate detector index (simplified)
            n = round(angleToTarget * f / delta_x_d);
            n = max(min(n, floor(Nd/2)), -floor(Nd/2));
            
            % Simulate active pixels (simplified)
            Active = Nd/2 - abs(n) + randn()*10;
            Active = max(min(Active, Nd), 0);
        end
        
        function updateHistoryBuffers(this, LOS, LOS_rate, Active, current_angle)
            % Update history buffers with new values
            this.LOS_hist = [LOS; this.LOS_hist(1:end-1)];
            this.LOS_rate_hist = [LOS_rate; this.LOS_rate_hist(1:end-1)];
            this.LOS_m_cmd_int_hist = [this.LOS_m_cmd_int; this.LOS_m_cmd_int_hist(1:end-1)];
            this.cmd_turn_rate_hist = [this.cmd_turn_rate; this.cmd_turn_rate_hist(1:end-1)];
            this.Active_hist = [Active/this.getNd(); this.Active_hist(1:end-1)];
            this.Active_rate_hist = [(Active - this.Active_old)/(this.getNd()*this.dt); this.Active_rate_hist(1:end-1)];
            this.current_angle_hist = [current_angle/(2*pi); this.current_angle_hist(1:end-1)];
        end
        
        function Pd = calculateDamageProbability(this)
            % Calculate damage probability
            St_prime = this.T_pos - this.S;
            Range = norm(St_prime);
            
            % Compute off-axis angle
            Vt = [-this.Vt_mag, 0];
            top = cross([this.V, 0], [Vt, 0]);
            sin_theta_off = top(3)/(norm(this.V)*norm(Vt));
            cos_psi = sin_theta_off;
            A = 2 * this.r_tgt_missile * this.L_target;
            Ae = A * cos_psi;
            if Ae < pi * this.r_tgt_missile^2
                Ae = pi * this.r_tgt_missile^2;
            end
            
            % Compute effective spreading area
            As = 2 * pi * (Range^2) * (1 - cos(this.phi));
            
            % Exponential fuze/damage probability
            Pd = 1 - exp(-this.K * (Ae/As));
        end
        
        function obs = getInitialObservation(this)
            % Get initial observation (all zeros)
            obs = zeros(this.inputSize, 1);
        end
        
        function obs = getLastObservation(this)
            % Get last observation (reuse previous values)
            obs = [
                %mean(this.LOS_hist);
                %mean(this.LOS_rate_hist);
                %mean(this.LOS_m_cmd_int_hist);
                mean(this.cmd_turn_rate_hist);
                mean(this.Active_hist);
                mean(this.Active_rate_hist);
                mean(this.current_angle_hist);
                this.Active_old/this.getNd();
                (this.Active_old - this.Active_hist(end))/(this.getNd()*this.dt);
                this.cmd_turn_rate/this.max_turn_rate
            ];
        end
        
        function Nd = getNd(this)
            % Get number of detectors
            Nd = 2.1998e+04;
        end
        
        function delta_x_d = getDeltaXd(this)
            % Get detector element width
            delta_x_d = 1.2500e-05;
        end
        
        function f = getFocalLength(this)
            % Get focal length
            f = 0.0030;
        end
    end
end