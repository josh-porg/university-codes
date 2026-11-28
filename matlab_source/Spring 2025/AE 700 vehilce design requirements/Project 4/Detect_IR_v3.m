% This function calculates the Detector Number and Number of Pixels illuminated
% in a 2-D Flyout Simulation

% The calling program must be intialized like this:
% V=[0 V_mag]; % Initial Interceptor Velocity; Format: [Vx Vy]
% S=[0 0]; % Initial Interceptor Position; Format: [Sx Sy]

% Programming Sequence
% 1. Interceptor Position and Velocity
% ---> 2. Interceptor Detection
% 3. Interceptor Turns
% 4. Interceptor Fuze, Pd?
% 5. Interceptor Moves
% 6. Target Moves

% Function requires INPUT of the following variables:
% V - Interceptor velocity vector [m/s]
% S - Interceptor position vector [m]
% wd - Width of the detector [m]
% f - focal length [m]
% Nd - Number of Detectors [~]
% Vt_mag - Target velocity magnitude [m/s]
% delta_x_d - width of one detector [m]
% time_elapsed - Time Elapsed by multiplying an increment (i=0,1,2...) by
% delta_t
% r_tgt_missile - Radius of target missile [m]
% L - Length of Motor Casing [m]
% x0_t - Target initial x-position at time_elapsed = 0 [+/- m]
% y0_t - Target initial y-position [+/- m]

% Function will OUTPUT the following variables:
% n - The integer number position of your detector
%   Note this is relative to center.
%   -1 is left first detector
%   +1 is right first detector from center
% If you have a center detector, it will report "0" for center
% Active - Number of illuminated detectors

function [n,Active] = Detect_IR_v3(V,S,Nd,wd,f,delta_x_d,time_elapsed,...
    r_tgt_missile,L,Vt_mag,x0_t,y0_t, T_pos)
% Target Initial Conditions
%St=[(-Vt_mag*time_elapsed)+x0_t y0_t]; % Initial Interceptor Position in meters
 St = T_pos;
St_prime = St - S; % Target Position Relative to Interceptor

% hold on
% % % plot(S(1),S(2),'ro')
% plot(St(1),St(2),'b*')

theta_need = acos(dot(V,St_prime)/(norm(V)*norm(St_prime))); % Needed Turn Angle
check = cross([St_prime 0], [V 0]); % Check of the Cross product sign indicating direction
if check(3) > 0
    theta_need = -theta_need;
else
end

theta_max = atan(wd/(2*f)); % Theta_max at the edge of the IR array. (Assuming corrected lens fall off)

if abs(theta_need)>= theta_max % Checks if the needed angle is beyond the detector width.
    n = 0;
    Active = 0;
    return
else
    n = (theta_need/theta_max)*(Nd/2);
    if theta_need < 0
        n = floor(n); % For negative theta_need the index for detector element is negative
    else
        n = floor(n) + 1; % For positive theta_need the index for detector element is positive
    end
end

Vt = [-Vt_mag 0];
top = cross([V 0],[Vt 0]);
sin_theta_off = top(3)/(norm(V)*norm(Vt)); % Generates projected area calculation. 
cos_psi=sin_theta_off; % Generates projected area calculation.
A = 2*r_tgt_missile*L; % Value of area normal to surface and defined as D * L of target missile.
Ae = A*cos_psi; % Projection of Area off normal to surface.
Le = L*cos_psi; % Projection of Length off normal to surface pixel illumination.

if Ae < pi*r_tgt_missile^2 % If the missile is in a "tail-chase"
    Le = 2*r_tgt_missile;
else
end

x_illum = (Le*f)/norm(St_prime); % Using the scale, or magnification, equation
Active = 2*(floor(x_illum/delta_x_d)+1)-1; % This assigns active pixel count.

return