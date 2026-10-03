% Parameters
Lx = 1.0;           % Domain length in x-direction (m)
Ly = 1.0;           % Domain length in y-direction (m)
Nx = 100;           % Number of grid points in x-direction
Ny = 100;           % Number of grid points in y-direction
dx = Lx / (Nx - 1); % Grid spacing in x-direction
dy = Ly / (Ny - 1); % Grid spacing in y-direction
dt = 0.001;         % Time step (s)
T = 1;            % Total simulation time (s)
S_L0 = 0.1;         % Unstretched laminar flame speed (m/s)
Markstein = 0.01;   % Markstein length (m)
D = 1e-5;           % Diffusion coefficient (m^2/s)

% Initialize grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Initial flame front (e.g., a circle)
R = 0.2;            % Initial radius of the flame front
phi = sqrt((X - Lx/2).^2 + (Y - Ly/2).^2) - R;

% Unsteady velocity field (e.g., sinusoidal perturbation)
u0 = 0.1;           % Base x-velocity (m/s)
v0 = 0.0;           % Base y-velocity (m/s)
omega = 2 * pi / T; % Frequency of unsteady flow

% Time-stepping loop
for t = 0:dt:T
    % Unsteady velocity field
    u = u0 * (1 + 0.1 * sin(omega * t) * Y / Ly);
    v = v0 * (1 + 0.1 * cos(omega * t) * X / Lx);
    
    % Compute gradients of phi
    [dphidx, dphidy] = gradient(phi, dx, dy);
    grad_phi = sqrt(dphidx.^2 + dphidy.^2);
    
    % Compute curvature (kappa)
    [d2phidx2, d2phidy2] = gradient(dphidx, dx, dy);
    kappa = (d2phidx2 + d2phidy2) ./ (grad_phi + 1e-6); % Avoid division by zero
    
    % Modify flame speed with curvature and stretch effects
    S_L = S_L0 * (1 - Markstein * kappa);
    
    % Update phi using the level set equation
    phi = phi - dt * (u .* dphidx + v .* dphidy) + dt * S_L .* grad_phi;
    
    % Diffusion term (optional, for smoothing)
    phi = phi + D * dt * del2(phi, dx, dy);
    
    % Visualization
    contourf(X, Y, phi, [0, 0], 'LineColor', 'r');
    title(['Time = ', num2str(t)]);
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal;
    drawnow;
end