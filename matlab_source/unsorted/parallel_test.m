%% parallel test
parpool();
disp(gcp)



tic
for i = 0:35
    factor(floor(9999999+exp(i+1)));
end
disp("factor time: "+toc);

tic
parfor i = 0:35
    factor(floor(9999999+exp(i+1)));
end
disp("factor time par: "+toc);

%% lorenz sim test 
Beta = [10; 28; 8/3]; % chaotic values
x0 = [0; 1; 20]; % initial condition
dt = 0.001; % times step
timeSpan = dt:dt:50; % time span

% set otions to alow for more accurate simulation and control of error tolerances
% set relative tolerance to 10^-12 and absoute tolerence of each component
% to the 10^-12
options = odeset('RelTol',1e-12, 'AbsTol',1e-12*ones(1,3));

% output is time and the time history of state x
% wrap lorenz function as function of x and t for ode45
% ode 45 wants time variable even though our DE is not time dependedent
tic
for i = 1:10
    [t,x] = ode45(@(t,x) lorenz(t,x,Beta), timeSpan, x0, options);
end
disp("lorenz time: "+toc);

tic
parfor i = 1:10
    [t,x] = ode45(@(t,x) lorenz(t,x,Beta), timeSpan, x0, options);
end
disp("lorenz time par: "+toc);

%%clean up
delete(gcp('nocreate'))