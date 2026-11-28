function IR_Seeker_Flyout_GUI()
    % --- Parameters ---
    V_mag = 400; % Interceptor speed (m/s)
    Vt_mag = 200; % Target speed (m/s)
    r_tgt_missile = 0.15; % Target missile radius (m)
    L = 0.6; % Target missile length (m)
    K = 90; % Number of fragments
    phi = deg2rad(30); % Fragment half-angle
    max_turn_rate = deg2rad(100); % rad/sec
    delta_t = 0.1;
    f = 0.1;
    Nd = 11;
    wd = 0.05;
    delta_x_d = wd / Nd;
    x0_t = 2000;
    y0_t = 0;

    % Initial Conditions
    S = [0; 0];
    V = [0; V_mag];
    positions = S';
    target_positions = [];
    Pd = 0;
    time_elapsed = 0;

    % --- Create Figure and UI ---
    fig = figure('Name','IR Seeker Flyout','NumberTitle','off');
    ax = axes('Parent',fig);
    hold on; grid on;
    axis equal;
    xlim([0 x0_t+100]);
    ylim([-500 500]);
    xlabel('X (m)'); ylabel('Y (m)');
    title('IR Seeker Flyout Animation');
    
    missile_line = plot(NaN, NaN, 'ro-', 'DisplayName','Interceptor');
    target_line = plot(NaN, NaN, 'b--', 'DisplayName','Target');
    missile_dot = plot(NaN, NaN, 'ro');
    target_dot = plot(NaN, NaN, 'bo');
    legend();

    % Pause Button
    uicontrol('Style', 'togglebutton', 'String', 'Pause', ...
        'Position', [20 20 100 30], 'Callback', @pause_callback, ...
        'Tag', 'pauseButton');

    function pause_callback(src, ~)
        if get(src, 'Value')
            set(src, 'String', 'Resume');
        else
            set(src, 'String', 'Pause');
        end
    end

    % --- Simulation Loop ---
    i = 1;
    while Pd < 0.9 && time_elapsed < 60
        time_elapsed = time_elapsed + delta_t;
        St = [x0_t - Vt_mag * time_elapsed; y0_t];
        St_prime = St - S;
        target_positions = [target_positions; St'];

        % IR Seeker Detection
        theta_need = acos(dot(V, St_prime) / (norm(V) * norm(St_prime)));
        cprod = cross([St_prime; 0], [V; 0]);
        if cprod(3) > 0
            theta_need = -theta_need;
        end

        theta_max = atan(wd / (2 * f));
        if abs(theta_need) > theta_max
            break; % No detection
        end

        % Turn Missile
        turn_angle = max(min(theta_need, max_turn_rate * delta_t), -max_turn_rate * delta_t);
        current_angle = atan2(V(2), V(1)) + turn_angle;
        V = V_mag * [cos(current_angle); sin(current_angle)];

        % Update Position
        S = S + V * delta_t;
        positions = [positions; S'];

        % Calculate Pd
        Vt = [-Vt_mag; 0];
        cprod2 = cross([V; 0], [Vt; 0]);
        sin_theta_off = cprod2(3) / (norm(V) * norm(Vt));
        cos_psi = sin_theta_off;
        A = 2 * r_tgt_missile * L;
        Ae = A * cos_psi;
        if Ae < pi * r_tgt_missile^2
            Ae = pi * r_tgt_missile^2;
        end
        As = 2 * pi * (norm(St - S)^2) * (1 - cos(phi));
        Pd = 1 - exp(-K * (Ae / As));

        % --- Animation Update ---
        set(missile_line, 'XData', positions(:,1), 'YData', positions(:,2));
        set(target_line, 'XData', target_positions(:,1), 'YData', target_positions(:,2));
        set(missile_dot, 'XData', S(1), 'YData', S(2));
        set(target_dot, 'XData', St(1), 'YData', St(2));
        drawnow;

        % Pause Logic
        while get(findobj('Tag','pauseButton'), 'Value')
            pause(0.1);  % Delay while paused
        end

        i = i + 1;
    end

    % Final Result
    disp(['Final Probability of Damage Pd = ', num2str(Pd)]);
end
