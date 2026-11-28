clear; close all; clc
% this is the true cessna 182
% note that insufficiently many lateral flap sections exist. at least 3
% spanwise flap sections are recommended per wing for desired lateral
% response.

dt = .01; % time step (seconds)
t_total = 20; % total simulation time (seconds)
angle = 5; % perturbed deflection (degrees)

disable_PAH = false; % use rigid flap?

k_command = .2; % commanded value of spring stiffness (N/m)
c_r_command = .05; % commanded value of rate damping ()
c_c_command = .005; % commanded value of coulomb damping ()
delt_0_command = deg2rad(5); % commanded value of deflection (rad)

I_flap = .9; % moment of inertia of flap about hinge line (kg*m^4)

x_0 = [k_command, c_r_command, c_c_command, delt_0_command, I_flap];

% best PAH values thus far:
x_0 = [73.3742, 5.5812, 0.6373, 0.1236, 28.2799]; % multiple times I have found this value

[y_all, time] = simulation(dt, t_total, angle, k_command, c_r_command, c_c_command, delt_0_command, I_flap, true);
plot_6dof_results(y_all, time, angle);

f = objective_function(x_0, dt, t_total, angle, true);

% [y_all, time] = simulation(dt, t_total, angle, k_command, c_r_command, c_c_command, delt_0_command, I_flap, false);
% plot_6dof_results(y_all, time, angle);
% f = objective_function(x_0, dt, t_total, angle, false);

% Define the objective function
obv_fun = @(x) objective_function(x, dt, t_total, angle, false);

% Define the initial design variables

% Define the bounds for the design variables
lb = [0, 0, .004, -2*pi, 1e-2];
ub = [1e2, 1e2, 1e2, 2*pi, 1e2];

% Simulated Annealing parameters
max_iter = 1000; % Number of iterations
initial_temp = 1; % Initial temperature
final_temp = 1e-3; % Final temperature
alpha = 0.99; % Cooling rate

% Initialize the solution
current_solution = lb + rand(1, numel(lb)) .* (ub - lb);
current_fitness = obv_fun(current_solution);
best_solution = current_solution;
best_fitness = current_fitness;

% Main SA loop
for iter = 1:max_iter
    temp = initial_temp * alpha^(iter - 1);
    
    % Generate a new candidate solution
    candidate_solution = current_solution + (rand(1, numel(lb)) - 0.5) .* (ub - lb) * 0.1;
    candidate_solution = max(min(candidate_solution, ub), lb);
    candidate_fitness = obv_fun(candidate_solution);
    
    % Acceptance probability
    delta = candidate_fitness - current_fitness;
    if delta < 0 || rand < exp(-delta / temp)
        current_solution = candidate_solution;
        current_fitness = candidate_fitness;
    end
    
    % Update the best solution found
    if current_fitness < best_fitness
        best_solution = current_solution;
        best_fitness = current_fitness;
    end
    
    % Display progress
    disp(['Iteration: ', num2str(iter), ', Current Fitness: ', num2str(current_fitness), ', Best Fitness: ', num2str(best_fitness)]);
    
    % Check if temperature is below the final temperature
    if temp < final_temp
        break;
    end
end

% Display the optimal solution
disp('Optimal design variables:');
disp(best_solution);
disp('Optimal objective function value:');
disp(best_fitness);

[y_all, time] = simulation(dt, t_total, angle, best_solution(1), best_solution(2), best_solution(3), best_solution(4), best_solution(5), false);
plot_6dof_results(y_all, time, angle);

function Xdot = aircraft(t, y, cv, stab_cont_derivs, I_flap, disablePAH)
    i = 1;
    states = y;
    controls = cv;
    
    VT      = states(0+i); %long
    alpha   = states(1+i); %long
    beta    = states(2+i);
    
    phi     = states(3+i);
    theta   = states(4+i); % long
    psi     = states(5+i);
    
    P       = states(6+i);
    Q       = states(7+i); %long
    R       = states(8+i);
    
    % Initial Control surfaces
    del_t_x   = states(12+i);
    del_e_x   = states(13+i);
    del_a_x   = states(14+i);
    del_r_x   = states(15+i);
    
    % PAH flap states
    % flap 1 (left)
    del_f_1  = states(16+i);
    del_f_dot_1  = states(17+i);
    % flap 2 (right)
    del_f_2  = states(18+i);
    del_f_dot_2  = states(19+i);
    
    % Control inputs
    del_t   = controls(0+i);
    del_e   = controls(1+i);
    del_a   = controls(2+i);
    del_r   = controls(3+i);
    
    % PAH control inputs
    k = controls(4+i);
    c_r = controls(5+i);
    delta_0 = controls(6+i);
    c_c = controls(7+i);
    
    u = VT*cos(alpha)*cos(beta);
    v = VT*sin(beta);
    w = VT*sin(alpha)*cos(beta);
    
    % Constants
    g = 32.14741; %ft/s2 %g*3.2808399;
    
    % Trim values
    alphatrim = 0; % convert input trim to rad
    Ptrim        = 0; %(rad/s)
    Qtrim        = 0; %(rad/s)
    Rtrim        = 0; %(rad/s)
    VTtrim       = 220.1; %(ft/s)
    utrim        = VTtrim*cos(alphatrim); %(ft/s)
    
    S    = 174; %(ft2)
    cbar = 4.9; %(ft)
    b    = 36; %(ft)
    AR   = 7.448;
    mass = 2650; %(lbm)
    I_xx = 948; %(slug-ft2)
    I_yy = 1346; %(slug-ft2)
    I_zz = 1967; %(slug-ft2)
    I_xz = 0.00; %(slug-ft2)
    
    % Flap parameters
    c_f = .25; % flap cord (m)
    C_h_0 = 0; % 0 angle of attack 0 deflection hinge moment (~)
    C_h_alpha = -.0584; % control surface hinge moment deriveative due to angle of attack (~)
    C_h_delta = -.363; % control surface moment derivative due to flap deflection (~)
    
    
    % location of Flap 1 (right wing)
    r_f_i_1 = [0; 0; 0]; % distance in body coordinates from the aircraft centre of gravity to the inboard most section of the control surface hingeline (m)
    r_f_o_1 = [0; b/2; 0]; % distance in body coordinates from the aircraft centre of gravity to the outboard most section of the control surface hingeline (m)
    r_f_i_TE_1 = [-c_f; 0; 0]; % distance in body coordinates from the aircraft centre of gravity to the trailing edge of the inboard most section of the undeflected control surface (m)
    r_f_cg_B_1 = [0; b/4; 0]; % distance in body coordinates from the aircraft centre of gravity to the location of the hingleline at the section of flap contining its centre of gravity (m)
    
    % location of flap 2
    r_f_i_2 = [0; -b/2; 0]; % distance in body coordinates from the aircraft centre of gravity to the inboard most section of the control surface hingeline (m)
    r_f_o_2 = [0; 0; 0]; % distance in body coordinates from the aircraft centre of gravity to the outboard most section of the control surface hingeline (m)
    r_f_i_TE_2 = [-c_f; -b/2; 0]; % distance in body coordinates from the aircraft centre of gravity to the trailing edge of the inboard most section of the undeflected control surface (m)
    r_f_cg_B_2 = [0; -b/4; 0]; % distance in body coordinates from the aircraft centre of gravity to the location of the hingleline at the section of flap contining its centre of gravity (m)
    
    
    
    % aircraft coefficients
    c_l_alpha_wf = 2*pi; % lift coefficient of wing fuselage combination per radian (1/radian)
    c_l_alpha_h = 2*pi; % lift coefficient of wing fuselage combination per radian (1/radian)
    tau_delta_f = .45; % efficiency of flap (~) from roskam stability and control part 1 page 61 figure 2.23
    eta_h = .9; % dynamic pressure ratio of horizontal tail to free stream (~)
    S_h = 1.75; % horizontal tail planform area (m^2)
    d_elpsilon_d_alpha_h = .111; % downwash gradient on the horizontal tail (1/radians) from AAA (can also compute using roskam part 6 page 262)
    %I_h = .9;
    I_h = I_flap;
    
    % Stop spring range and coeffcient
    delta_min = -deg2rad(20); % minimum deflection angle (radians)
    delta_max = deg2rad(20); % maximum deflection angle (radians)
    k_stop = 1000; % deflection limiting hardening spring stiffness N/(m^2)
    
    % Note to self: flip flap limits for left aeileron because positive and
    % negative may be flipped
    
    % steady state
    
    C_L1  = 0.307;
    C_m1  = 0.0;
    C_mt1 = 0.0;
    
    % S&C derivatives
    C_D0bar   = 0.032;
    e=.75;
    
    C_yda     = 0;
    C_Lq      = 3.9;
    C_Ladot   = 1.7;
    C_mtu     = 0;
    C_mtalpha = 0;
    C_yb      = -.393;      C_yp      = -.075;
    C_yr      = .214;       C_La      = 4.41;
    C_Lu      = 0;          C_lb      = -.0923;
    C_lp      = -.484;      C_lr      = .0798;
    C_ma      = -.613;      C_mq      = -12.4;
    C_madot   =-7.27 ;      C_mu      = 0;
    C_nb      = .0587;      C_np      = -.0278;
    C_nr      = -.0937;     C_ydr     = .187;
    C_ldela   = .229;       C_ldelr   = .0147;
    C_ndela   = -.0216;     C_ndelr   = -.0645;
    C_Ldele   = .43;        C_mdele   = -1.122;
    C_yda     = 0;
    C_mtu     = 0;
    C_mtalpha = 0;
    
    mass = mass/g; % lbm to slug
    rho  = 0.002048;
    qbar = 0.5*rho*(VT^2);
    
    % aero forces and moment
    C_L  = C_L1 + C_La*(alpha - alphatrim) + (C_Lq*(Q - Qtrim)*(cbar/2))/utrim + (C_Lu*(u - utrim))/utrim + C_Ldele*(del_e_x ) + (1/2) * C_La* tau_delta_f * del_f_1 + (1/2) * C_La* tau_delta_f * del_f_2 ; % be careful with sign of flap deflection for left flap
    C_D  = C_D0bar + (C_L^2)/(pi*AR*e); % flap drag not implemented
    
    C_Y  = C_yb*beta + (C_yp*(P - Ptrim)*(b/2))/utrim + (C_yr*(R - Rtrim)*(b/2))/utrim + C_yda*del_a_x + C_ydr*del_r_x; % flap induced adverse yaw not modeled
    C_ls = C_lb*beta + (C_lp*(P - Ptrim)*(b/2))/utrim + (C_lr*(R - Rtrim)*(b/2))/utrim + C_ldela*del_a_x + C_ldelr*del_r_x - (1/2) * C_La* tau_delta_f * del_f_1 * r_f_cg_B_1(2) - (1/2) * C_La* tau_delta_f * del_f_2 * r_f_cg_B_2(2);
    % note: flap moment not implemented - considered negigible
    C_ms = C_m1 + C_ma*(alpha - alphatrim) + (C_mq*(Q - Qtrim)*(cbar/2))/utrim + (C_mu*(u - utrim))/utrim + C_mdele*(del_e_x ) + ...
        2*(C_m1*(u - utrim))/utrim + (C_mtu + 2*C_mt1)*(u - utrim)/utrim + C_mtalpha*(alpha - alphatrim);
    % flap yaw moments not implemented
    C_ns = C_nb*beta + (C_np*(P - Ptrim)*(b/2))/utrim + (C_nr*(R - Rtrim)*(b/2))/utrim + C_ndela*(del_a_x) + C_ndelr*del_r_x;
    C_xa = C_L*sin(alpha) - C_D*cos(alpha);
    C_ya = C_Y;
    C_za = -C_L*cos(alpha) - C_D*sin(alpha);
    C_l  = C_ls*cos(alpha) - C_ns*sin(alpha);
    C_m  = C_ms;
    C_n  = C_ls*sin(alpha) + C_ns*cos(alpha);
    
    %Thrust model
    %???????????????????????????????????????????????????????????????????????????????????????????????
    T = 322.5669;
    
    % grav accel
    g_x = -g*sin(theta);
    g_y = g*sin(phi)*cos(theta);
    g_z = g*cos(phi)*cos(theta);
    
    % body accel
    a_x = (qbar*S*C_xa + T)/mass;
    a_y = (qbar*S*C_ya)/mass;
    a_z = (qbar*S*C_za)/mass;
    
    % moment
    L = C_l*qbar*S*b;
    M = C_m*qbar*S*cbar;
    N = C_n*qbar*S*b;
    % Differential equations:
    % Body linear velocities
    udot = a_x + g_x + R*v - Q*w;
    vdot = a_y + g_y - R*u + P*w;
    wdot = a_z + g_z + Q*u - P*v;
    
    % Total velocity and airflow angles
    VTdot    = (u*udot + v*vdot + w*wdot)/VT;
    alphadot = (u*wdot - w*udot)/((u^2) + (w^2));
    betadot  = ((VT*vdot - v*VTdot)/((u^2) + (w^2)))*cos(beta);
    %betadot  = (vdot*((u^2) + (w^2)) - v*(u*udot + w*wdot))/((VT^2)*(sqrt((u^2) + (w^2))));
    
    
    
    % Body angular velocities
    Pdot = (I_zz*L + I_xz*N - (I_xz*(I_yy - I_xx - I_zz)*P + ((I_xz^2) + I_zz*(I_zz - I_yy))*R)*Q)/(I_xx*I_zz - (I_xz^2));
    Qdot = (M - (I_xx - I_zz)*P*R - I_xz*((P^2) - (R^2)))/I_yy;
    Rdot = (I_xz*L + I_xx*N + (I_xz*(I_yy - I_xx - I_zz)*R + ((I_xz^2) + I_xx*(I_xx - I_yy))*P)*Q)/(I_xx*I_zz - (I_xz^2));
    
    % Inertial (Euler) angles
    phidot   = P + (R*cos(phi) + Q*sin(phi))*tan(theta);
    thetadot = Q*cos(phi) - R*sin(phi);
    psidot   = (Q*sin(phi) + R*cos(phi))/cos(theta);
    
    % Inertial position dynamics
    x_Idot = (cos(theta)*cos(psi))*u + (-cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi))*v + (sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi))*w;
    y_Idot = (cos(theta)*sin(psi))*u + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*v + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*w;
    z_Idot = (-sin(theta))*u + (sin(phi)*cos(theta))*v + (cos(phi)*cos(theta))*w;
    
    % Control surface dynamics
    Athrottle=-7;
    Bthrottle=7;
    Aelevator=-30;
    Belevator=30;
    Aaileron=-30;
    Baileron=30;
    Arudder=-30;
    Brudder=30;
    
    del_tdot = Athrottle*del_t_x + Bthrottle*del_t;
    del_edot = Aelevator*del_e_x + Belevator*del_e;
    del_adot = Aaileron*del_a_x + Baileron*del_a;
    del_rdot = Arudder*del_r_x + Brudder*del_r;
    
    % get flap angle
    q_bar = 1/2 * rho * VT^2; % dynamic pressure (Pa)
    % 
    % determine prefered damping constant
    % set c_r for airfoil responce properties
    zeta = .7; % damping ratio
    mu = 4*I_h* (q_bar * c_f^2 * C_h_delta - k); % airfoil responce parameter
    % set damping to make that true (note culombic damping will be suboptimal
    c_d = 2*I_h*zeta * sqrt(-mu * (2*I_h*zeta - 1) * (2*I_h*zeta + 1)) / (4*I_h^2*zeta^2-1); % total damping coefficient
    % unsure if i should use the one that acounts for culombic damping (good
    % near fp) or the one that doesn't (good far from f.p.) I think that as
    % this is transient responce culombic damping should be acounted for
    %c_r = max(c_d - c_c,0); % for now I will comment out this to disable 
    
    
    % Adaptive_Airfoil_Dynamics is outdated 
    %df = Adaptive_Airfoil_Dynamics(t, [del_f, del_f_dot], I_h, k, delta_0, c_r, c_c, c_f, C_h_0, C_h_alpha, C_h_delta, alpha, q_bar);
    df_1 = Adaptive_Aerocompliant_Airfoil_Dynamics(t, [del_f_1, del_f_dot_1], I_h, k, delta_0, c_r, c_c, c_f, C_h_0, C_h_alpha, C_h_delta, delta_min, delta_max, k_stop, r_f_i_1, r_f_o_1, r_f_i_TE_1, r_f_cg_B_1, rho, [u;v;w], [P;Q;R]);
    del_f_dot_1 = df_1(1);
    del_f_dotDot_1 = df_1(2);
    
    df_2 = Adaptive_Aerocompliant_Airfoil_Dynamics(t, [del_f_2, del_f_dot_2], I_h, k, delta_0, c_r, c_c, c_f, C_h_0, C_h_alpha, C_h_delta, delta_min, delta_max, k_stop, r_f_i_2, r_f_o_2, r_f_i_TE_2, r_f_cg_B_2, rho, [u;v;w], [P;Q;R]);
    del_f_dot_2 = df_2(1);
    del_f_dotDot_2 = df_2(2);
    
    
    if disablePAH
        del_f_dot_1 = 0;
        del_f_dot_2 = 0;
    end
    
    Xdot = [VTdot; alphadot; betadot; phidot; thetadot; psidot; Pdot; ...
        Qdot; Rdot; x_Idot; y_Idot; z_Idot; del_tdot; del_edot; del_adot; del_rdot; del_f_dot_1; del_f_dotDot_1; del_f_dot_2; del_f_dotDot_2];
end

function [y_all, time] = simulation(dt, t_total, angle, k_command, c_r_command, c_c_command, delt_0_command, I_flap, disable_PAH)
    % define IC and controls
    time = [0:dt:t_total]';
    t_0 = 1;
    t_end = length(time);
    
    throttle = zeros(size(time));
    elevator = zeros(size(time));
    aileron = zeros(size(time));
    rudder = zeros(size(time));
    elevator(1/dt:1.5/dt) = ones(size(1/dt:1.5/dt)).*angle.*pi./180; %singlet
    %elevator(1.5/dt:2/dt) = -ones(size(1.5/dt:2/dt)).*angle.*pi./180; % add to make doublet
    aileron(1/dt:1.5/dt) = ones(size(1/dt:1.5/dt)).*angle.*pi./180;


    % stockastic perturbations

    % Generate random noise
    n_samples = length(time);
    noise = randn(n_samples,1);

    % Apply a moving average filter to smooth the noise
    window_size = 50; % Size of the moving window
    smooth_noise = filter(ones(1, window_size) / window_size, 1, noise);
    % filiter again on the fine scale to make it more differentiable
    window_size = 5; % Size of the moving window
    smooth_noise = filter(ones(1, window_size) / window_size, 1, smooth_noise);
    
    %apply to surface
    %elevator = smooth_noise;
    %aileron = smooth_noise;


    % PAH commands
    k = ones(size(time)) .* k_command; % PAH spring const
    c_r = ones(size(time)) .* c_r_command; % rate damping coefficient (N*s/rad^2)
    c_c = ones(size(time)) .* c_c_command; % coulomb damping coefficient (N*s/rad^2)
    delta_0 = ones(size(time)) .* delt_0_command; % commanded deflection (radians)
    % colombic damping coefficient (N)pm
    
    %k = k + .03*randn(size(time));
    %delta_0 = delta_0 + .003*randn(size(time));
    
    
    % elevator(1.5/dt:2/dt) = ones(size(1.5/dt:2/dt)).*angle.*-pi./180.;
    % elevator = elevator  + .003*randn(size(elevator));
    VT = zeros(size(time));
    VT(1) = 220.1;
    Alphaest = zeros(size(time));
    Betaest = zeros(size(time));
    phi = zeros(size(time));
    theta = zeros(size(time));
    psi = zeros(size(time));
    P = zeros(size(time));
    Q = zeros(size(time));
    R = zeros(size(time));
    ned = [zeros(size(time)),zeros(size(time)),zeros(size(time))];
    delta_f_1 = 0; % inital flap 1 position
    delta_f_dot_1 = 0; % initial flap 1 rate
    delta_f_2 = 0; % inital flap 2 position
    delta_f_dot_2 = 0; % initial flap 2 rate
    
    cv = [throttle,   elevator ,   aileron  ,   rudder, k, c_r, delta_0, c_c];
    X_0 = [VT(1), Alphaest(1), Betaest(1), phi(1), theta(1), psi(1), ...
        P(1), Q(1), R(1), ned(1,1), ned(1,2), -ned(1,3), throttle(t_0), elevator(t_0) , ...
        aileron(t_0) , rudder(t_0),...
        delta_f_1, delta_f_dot_1, delta_f_2, delta_f_dot_2];
    
    y = X_0;
    
    %initialize y_all and t_all with IC as first element
    y_all = zeros([length(t_0:t_end),20]);       t_all = zeros([length(t_0:t_end),1]);
    y_all(1,:) = y;     t_all(1) = time(t_0);       J = 2;
    
    for i = (t_0):(length(time)-1)
        cv1 = cv(i,:); %update control inputs and time span for ODE45 each iteration
    
        % PID controller for the PAH
        %k(demanded)
    
        tspan = [time(i) time(i+1)];       dt = time(i+1) - time(i);
        %disp("current time t = " + i*dt)
        %[t,y] = ode45(@(t,y) aircraft(t, y, cv1, [], I_flap, disable_PAH), tspan, X_0);
        [t,y] = RK4(@(t,y) aircraft(t, y, cv1, [], I_flap, disable_PAH), tspan, dt, X_0');
        y = y';

        X_0 = y(size(y, 1),:);  % set new IC as ending state
        t_all(J) = t(end);       y_all(J,:) = y(end,:);
        J = J + 1; % indexer for recording y and t
    end
end

function f = objective_function(x, dt, t_total, angle, disable_PAH)
    % unpack x
    k_command = x(1); c_r_command = x(2); c_c_command = x(3); delt_0_command = x(4); I_flap = x(5);

    num_samples = 1; % use multiple for stochasitc or 1 for detministic pertubations

    for i = 1:length(num_samples)
        % run simulation
        [y_all, time] = simulation(dt, t_total, angle, k_command, c_r_command, c_c_command, delt_0_command, I_flap, disable_PAH);
    
        % i dont think i need to do this but ben does it so i will keep it
        %13-16 are elevator alteroion etc.
        var_guid = {'U (ft/s)', '$\alpha (^\circ)$', '$\beta (^\circ)$', '$\phi (^\circ)$' '$\theta (^\circ)$', '$\psi (^\circ)$', 'P ($^\circ / s$)', 'Q ($^\circ / s$)','R ($^\circ / s$)','N' ,'E' ,'D', 'throttle', 'elevator', 'aileron', 'rudder', '$\delta_{f,1}$', '$\dot{\delta_{f_1}}$','$\delta_{f_2}$', '$\dot{\delta_{f_2}}$'};
        for ii = 1:length(var_guid)
            Vname{ii} = y_all(:,ii);
        end
    
        % we will minimize maximum total acceleration
        f_accel_max = max(sqrt( diff(diff(Vname{10})).^2 + diff(diff(Vname{11})).^2 + diff(diff(Vname{12})).^2 ));
        f_vel = sum(sqrt(diff(diff(Vname{11})).^2 + diff(diff(Vname{12})).^2 ));

        %TODO: minimi
    
        f_accel_sum = sum(sqrt( diff(diff(Vname{10})).^2 + diff(diff(Vname{11})).^2 + diff(diff(Vname{12})).^2 ));
        f_accel_ave  = f_accel_sum / (length(Vname{10})-2);

        f(i) = f_accel_ave;
        f(i) = f_accel_max;
    end


    f = sum(f)/length(f);

    disp("f " + f)
    %plot_6dof_results(y_all, time, angle)
end

function plot_6dof_results(y_all, time, angle)
    %Convert to deg
    y_all(:,2:end) = y_all(:,2:end)*180/pi;
    
    %13-16 are elevator alteroion etc.
    var_guid = {'U (ft/s)', '$\alpha (^\circ)$', '$\beta (^\circ)$', '$\phi (^\circ)$' '$\theta (^\circ)$', '$\psi (^\circ)$', 'P ($^\circ / s$)', 'Q ($^\circ / s$)','R ($^\circ / s$)','N' ,'E' ,'D', 'throttle', 'elevator', 'aileron', 'rudder', '$\delta_{f,1}$', '$\dot{\delta_{f_1}}$','$\delta_{f_2}$', '$\dot{\delta_{f_2}}$'};
    
    % Plot Longitudinal variables (VT, VT_guid, theta, theta_guid, Q, del_t, del_e)
    for ii = 1:length(var_guid)
        Vname{ii} = y_all(:,ii);
    end
    
    figure
    ii = 0;
    for i = [2,8,12,17] %7,9
        ii = ii + 1;
        subplot(4,1,ii)
        plot(time,Vname{i}), grid on, title(var_guid{i}, 'Interpreter','latex'), xlabel('time', 'Interpreter','latex')
        %title(['Doublet Input Responses for $\delta_x$ = $\pm$', num2str(angle), '$^{\circ}$'], 'Interpreter', 'Latex')
    end
    
    figure
    plot(time(1:end-2),diff(diff(Vname{12}))), grid on, title("d^2"+var_guid{12}+"/dt^2", 'Interpreter','latex'), xlabel('time', 'Interpreter','latex')
    title(['Acceleration Doublet Input Responses for $\delta_x$ = $\pm$', num2str(angle), '$^{\circ}$'], 'Interpreter', 'Latex')
    
    % latteral plots
    figure
    ii=0;
    for i = [3 4 6 7 9, 15]
        ii = ii + 1;
        subplot(3,2,ii)
        plot(time,Vname{i}), grid on, title(var_guid{i}, 'Interpreter','latex'), xlabel('time', 'Interpreter','latex')
    end
    sgtitle(['Lat. Doublet Input Responses for $\delta_a$ = ', num2str(angle), '$^{\circ}$'], 'Interpreter', 'Latex')
    
    
    % PAH plots
    figure
    ii=0;
    for i = [17,18,19,20]
        ii = ii + 1;
        subplot(3,2,ii)
        plot(time,Vname{i}), grid on, title(var_guid{i}, 'Interpreter','latex'), xlabel('time', 'Interpreter','latex')
    end
end
