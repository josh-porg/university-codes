% This function calculates a detection presence and Range
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
% diameter - Diameter of the radar aperture [m]
% time_elapsed - Time Elapsed by multiplying an increment (i=0,1,2...) by
% lambda - Wavelength of the radar [m]
% x0_t - Target initial x-position [+/- m]
% y0_t - Target initial y-position [+/- m]
% Vt_mag - Target velocity magnitude [m/s]

% Function will OUTPUT the following variables:
% n - "0" if no return, "1" if return
% Range - If no return, sets Range = 'No Return', otherwise a value [m]

function [n,Range] = Detect_Radar(V,S,diameter,time_elapsed,lambda,x0_t,y0_t,Vt_mag)

% Target Initial Conditions
St=[(-Vt_mag*time_elapsed)+x0_t y0_t]; % Initial Interceptor Position in meters
St_prime = St - S; % Target Position Relative to Interceptor

% hold on
% % % plot(S(1),S(2),'ro')
% plot(St(1),St(2),'b*')

theta_need = acos(dot(V,St_prime)/(norm(V)*norm(St_prime))); % Needed Turn Angle
check = cross([St_prime 0], [V 0]);
if check(3) > 0
    theta_need = -theta_need;
else
end

theta_max = 1.3*lambda/diameter; % Theta_max at the edge of the beamwidth. (Assuming Good Power)

if abs(theta_need)>= theta_max % Checks if the needed angle is beyond the detector width.
    n = 0;
    Range = 'No Return';
    return
else
    n = 1;
    Range = norm(St_prime);
end

return