function [t,x] = getBlackboxData2D(x10,x20,x30)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Beta = [10; 28; 8/3]; % chaotic values
x0 = [456; 100; 29]; % initial condition
dt = 0.001; % times step
timeSpan = dt:dt:50; % time span

% set otions to alow for more accurate simulation and control of error tolerances
% set relative tolerance to 10^-12 and absoute tolerence of each component
% to the 10^-12
options = odeset('RelTol',1e-12, 'AbsTol',1e-12*ones(1,3));

% output is time and the time history of state x
% wrap lorenz function as function of x and t for ode45
% ode 45 wants time variable even though our DE is not time dependedent
[t,x] = ode45(@(t,x) blackBox2(t,x,Beta), timeSpan, x0, options);
end



