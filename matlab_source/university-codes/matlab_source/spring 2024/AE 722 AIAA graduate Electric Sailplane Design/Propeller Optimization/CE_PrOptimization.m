clear all
close all


V_flt_kts = [40:5:80];
%V_flt_kts = linspace(40,80, 50)

r_max_gear = 2.05105/2;


for i = 1:length(V_flt_kts)

    %radius1 = [400:100:1600]*10^(-3);
    radius1 = [1000:200:1200]*10^(-3);
    %radius1 = [700:100:1300]*10^(-3);
    %radius1 = [linspace(700e-3,r_max_gear,24), linspace(r_max_gear+1e-3,1300e-3,27)]; 

    for j = 1:length(radius1)
        Eta_Final(i,j) = Multi_Radius_Opt(V_flt_kts(i), radius1(j));
    end
end


figure
for i = 1:length(V_flt_kts)
    
    plot(radius1,Eta_Final(i,:), 's-'), hold on

end 

% save the variavles




save("Optimization_Prop_and_Velocity_13", "Eta_Final","radius1","V_flt_kts");

function Eta_Final = Multi_Radius_Opt(V_flt_kts, radius)

Theta_0 = deg2rad(45); % root incidence angle (radians)
Theta_1 = deg2rad(15); % ? (radians)
Theta_tw = deg2rad(-42); % ? (radians)
Theta_r_twistExponent = -0.7; % ? (~)
c_r = 0.1; % root cord (m)
P_desired = 14250; % W
P_desired = 1.9700e+04; % W
A_mean = [Theta_0, Theta_1, Theta_tw, Theta_r_twistExponent, c_r, P_desired];
A_old = A_mean;
eta_last = .1; 

% Set initial standard deviation for A
MEAN_std_devA = [0.02, .01, .01, 0.01, .02, 0]*3;

std_devA = MEAN_std_devA;

% Set number of variants (simulation sample size / number of particles(derivative combinations))
num_variants = 300; % use more for the first iteration
num_main = 300; % primary number of variants

% Initialize storage matrix for derivatives of each variant and PQR
eta_all = zeros(num_variants, 1);
Omega_all = zeros(num_variants, 1);
alphaL_all = zeros(num_variants, 100);

mse = zeros(num_variants, 1);
min_mse_all = 10; % make a matrix for the min mse of each iteration (cost tracking)
jjj = 0;

num_runs = 18;
% Initialize arrays to store variables for plotting
Theta_0_array = zeros(num_runs, 1);
Theta_1_array = zeros(num_runs, 1);
Theta_tw_array = zeros(num_runs, 1);
Theta_r_twistExponent_array = zeros(num_runs, 1);
c_r_array = zeros(num_runs, 1);

%figure

for iiii = 1:num_runs
    jjj = jjj + 1;

    parfor j = 1:num_variants
        % Apply random normal error to A
        A = Normal_Dist_AB(A_mean, std_devA);
        [eta_all(j), ~, ~, Omega_all(j), alphaL] = PropDesign_SetPower(A, V_flt_kts, radius);
        A_all(j, :) = A;
        alphaL_all(j, :) = alphaL;
        if any(alphaL(40:end)*180/pi > 16.9)
            eta_all(j) = 0;
        end
    end
    A_all = [A_all; A_old]; % add previous best

    % Find candidates: best fit lines (using mse?? top ~10% ??)
    [~, sorted_indx] = sort(eta_all, 'descend');
    cand_index = sorted_indx(1:round(0.60 * length(mse)));


    k_A = 1;
    A_old = A_mean;

    A_mean = A_all(cand_index(1), :);
    alphaL_all = alphaL_all(cand_index(1), :);
    std_devA = abs(std(A_all(cand_index, :), 0, 1) ./ A_old) .* k_A; % Specify dimension 3 for z direction
    A_old = A_all(cand_index, :);


    Theta_0 = A_mean(1) * 180 / pi;
    Theta_1 = A_mean(2) * 180 / pi;
    Theta_tw = A_mean(3) * 180 / pi;
    Theta_r_twistExponent = A_mean(4);
    c_r = A_mean(5);



    plot(iiii, max(eta_all)*100, 'b-o'), hold on
    title('Eta')
    xlabel('Iteration')
    ylabel('%')

    drawnow
    
    % Store variables for plotting
    Theta_0_array(iiii) = Theta_0;
    Theta_1_array(iiii) = Theta_1;
    Theta_tw_array(iiii) = Theta_tw;
    Theta_r_twistExponent_array(iiii) = Theta_r_twistExponent;
    c_r_array(iiii) = c_r;

    if (max(eta_all)-eta_last)/eta_last < .2e-6
        break
    end
    eta_last = max(eta_all);
end

% Plotting the variables
% figure
% subplot(3, 2, 1)
% plot( Theta_0_array)
% title('Theta_0')
% xlabel('Iteration')
% ylabel('Value')
% 
% subplot(3, 2, 2)
% plot(Theta_1_array)
% title('Theta_1')
% xlabel('Iteration')
% ylabel('Value')
% 
% subplot(3, 2, 3)
% plot( Theta_tw_array)
% title('Theta_tw')
% xlabel('Iteration')
% ylabel('Value')
% 
% subplot(3, 2, 4)
% plot(Theta_r_twistExponent_array)
% title('Theta_r_twistExponent')
% xlabel('Iteration')
% ylabel('Value')
% 
% subplot(3, 2, 5)
% plot(c_r_array)
% title('c_r')
% xlabel('Iteration')
% ylabel('Value')
% 
% subplot(3, 2, 6)
% plot(alphaL_all)
% title('Limited alpha')
% xlabel('Iteration')
% ylabel('Value')

Theta_0 = A_mean(1)*180/pi;
Theta_1 = A_mean(2)*180/pi;
Theta_tw = A_mean(3)*180/pi;
Theta_r_twistExponent = A_mean(4);
c_r = A_mean(5);
Omega = Omega_all(end) * 60 / (2*pi);
Eta_Final = max(eta_all);

% Print variables in a table-like format
fprintf('Variable Name\t\tValue\n');
fprintf('--------------------------------\n');
fprintf('Theta_0\t\t\t%.2f\n', Theta_0);
fprintf('Theta_1\t\t\t%.2f\n', Theta_1);
fprintf('Theta_tw\t\t\t%.2f\n', Theta_tw);
fprintf('Theta_r_twistExponent\t%.2f\n', Theta_r_twistExponent);
fprintf('c_r\t\t\t\t%.2f\n', c_r);
fprintf('Omega\t\t\t\t%.2f\n', Omega);
fprintf('Eta_Final\t\t\t%.2f\n', Eta_Final);

end

function A_or_B_norm = Normal_Dist_AB(A_or_B, std_dev)
    % st_dev is to be a matrix the same size as A or B with the std dev
    % desired for each term in that term's respective position

    % Apply normal distribution to each term in A_or_B, bounded within given limits
    A_or_B_norm = A_or_B  +  std_dev .* A_or_B  .* randn(size(A_or_B));

    % Define lower and upper bounds for each variable
    %  [Theta_0, Theta_1, Theta_tw, Theta_r_twistExponent, c_r, P_desired]
    lower_bounds = [20*pi/180, 8*pi/180, -65*pi/180, -1.5, .065, -Inf]; % Define lower bounds for each variable
    upper_bounds = [65*pi/180, 30*pi/180, -30*pi/180, -.6, 0.15, Inf]; % Define upper bounds for each variable

    % Apply bounds to each term in A_or_B_norm
    for i = 1:numel(A_or_B_norm)
        if A_or_B_norm(i) < lower_bounds(i)
            A_or_B_norm(i) = lower_bounds(i);
        elseif A_or_B_norm(i) > upper_bounds(i)
            A_or_B_norm(i) = upper_bounds(i);
        end
    end
end
