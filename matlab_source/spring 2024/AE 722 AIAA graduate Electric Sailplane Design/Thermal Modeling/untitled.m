clear all; close all; clc;

e = .9;
r = 0:175; % Turn radius m
A = 22.82;%linspace(10,60,length(r)); % Aspect Ratio
rho = 1.0834; % Density at a+*lt
%[~, ~, ~, rho] = atmosisa(convlength(3000, 'ft', 'm'));

%W2S = 450;%linspace(250:550,length(r)); % Wing Loading N/m^2
W2S = 280; % to mathc figure 89 of Thomas
%W2S = 450 % messing with wing loading


L_Dmax = 22.82;

C_L_max = (pi * A * e) / (2 * L_Dmax); % Maximum lift coefficient
%C_L_max = 2.7

C_D_max = 2 * (C_L_max^2)/(pi * A * e); % Maximum drag coefficient
C_D_0 = C_D_max / 2; % Parasite drag coefficient

%for our aircraft
C_D_0 = .0310;
C_L_max = 2.7;

n_PlotPoints = 100;

W2S_min_plot = 100;
W2S_max_plot = 700;
A_min_plot = 10;
A_max_plot = 40;

W2Ss = linspace(W2S_min_plot, W2S_max_plot, n_PlotPoints)
As = linspace(A_min_plot, A_max_plot, n_PlotPoints)

for i = 1:length(W2Ss)
    disp(i)
    W2S = W2Ss(i);
    V_ave_quast(i,:) = arrayfun(@(A) AverageCrossCountrySpeed_Quast(e, A, C_D_0, W2S, rho, C_L_max), As);
end

[W2Ssgrid, Asgrdid] = meshgrid(W2Ss,As);
%contour3(W2Ssgrid, Asgrdid, V_ave_quast, 1000, "LineWidth",.001), colormap("winter"); set(gca, "Color", [1,1,1]); colorbar;
contourf(W2Ssgrid, Asgrdid, V_ave_quast, 100, "LineWidth",.001), colormap("winter"); set(gca, "Color", [1,1,1]); colorbar;
title("average cross country speed as a function of wingloadina s aspect ratio");
xlabel("$Wing Loading, \frac{W}{S} (\frac{N}{m^2})$",Interpreter="latex");
zlabel("Average Cross Country Speed, V_ave, (m/s)",Interpreter="latex");
ylabel("Aspect Ratio, A, (\sim)",Interpreter="latex");

% for i = 1:length(W2Ss)
%     W2S = W2Ss(i);
%     V_ave_quast(i) = AverageCrossCountrySpeed_Quast(e, A, C_D_0, W2S, rho, C_L_max);
% end
% plot(W2Ss, V_ave_quast); title("The optimum wingloading to maximise cross country speed",Interpreter="tex"); xlabel("$Wing Loading, \frac{W}{S} (\frac{N}{m^2})$",Interpreter="latex"), ylabel("Average Cross Country Speed, V_ave, (m/s)")
