clear
% good typical results
% load flightData_Rigid_v0.mat
% [rms_total_control, vdv_total_control, c_f_control] = ride_quality_analysis(flight_data, true);
% load flightData_AC_v0.mat
% [rms_total, vdv_total, c_f] = ride_quality_analysis(flight_data, true);


% load flightData_Rigid.mat
% [rms_total_control, vdv_total_control, c_f_control] = ride_quality_analysis(flight_data, true);
% 
% load flightData_AC.mat
% [rms_total, vdv_total, c_f] = ride_quality_analysis(flight_data, true);

load Dryden_opt_flightData_Rigid.mat
load flightData_Rigid_v0.mat
load flightData_Rigid.mat
[rms_total_control, vdv_total_control, c_f_control] = ride_quality_analysis(flight_data, true);
load Dryden_opt_flightData_AC.mat
load flightData_AC_v0.mat
load flightData_UnOpt.mat
load flightData_AC.mat
[rms_total, vdv_total, c_f] = ride_quality_analysis(flight_data, true);


% load flightData_UnOpt.mat
% [rms_total, vdv_total, c_f] = ride_quality_analysis(flight_data, true);

vdv_total_control
vdv_total
FOM_control = rms_total_control + vdv_total_control + sum(c_f_control)
FOM = rms_total + vdv_total + sum(c_f)

%now superimpose the plots 
figure
load flightData_Rigid.mat
load Dryden_opt_flightData_Rigid.mat
load flightData_Rigid_v0.mat
load flightData_Rigid.mat
% compute the ride quality of hte 
ax = flight_data.ax(:); % Ensure column vector
ay = flight_data.ay(:);
az = flight_data.az(:);
fs = flight_data.fs;
t = flight_data.t;


% Apply frequency weighting with proper filter syntax
% Define numerator & denominator for Wk
numWk = [1 32 1859 793.9 81.89];
denWk = [1 2261 7007 20640 79.96];
% Analog to discrete at sample rate fs
fs = flight_data.fs;  % e.g., 200 Hz
Ts = 1/fs;
[bdWk, adWk] = bilinear(numWk, denWk, fs);

% Define numerator & denominator for Wd
numWd = [0 14.55 6.026 7.714];     % leading zero for proper order
denWd = [1 51.65 47.66 15.02];
[bdWd, adWd] = bilinear(numWd, denWd, fs);

% Filter your acceleration signals
az_weighted = filter(bdWk, adWk, az);     % vertical axis
ax_weighted = filter(bdWd, adWd, ax);     % longitudinal
ay_weighted = filter(bdWd, adWd, ay);     % lateral

% Remove initial transient (1 second)
transient = round(1*fs);
t = t(transient:end);
ax_weighted = ax_weighted(transient:end);
ay_weighted = ay_weighted(transient:end);
az_weighted = az_weighted(transient:end);

figure
% Calculate PSD using Welch's method
nfft = 2^nextpow2(length(t)/10);
window = hann(nfft);
noverlap = nfft/2;
[pxx_x, f] = pwelch(ax, window, noverlap, nfft, fs);
[pxx_xw, ~] = pwelch(ax_weighted, window, noverlap, nfft, fs);

[pxx_y, ~] = pwelch(ay, window, noverlap, nfft, fs);
[pxx_yw, ~] = pwelch(ay_weighted, window, noverlap, nfft, fs);

[pxx_z, ~] = pwelch(az, window, noverlap, nfft, fs);
[pxx_zw, ~] = pwelch(az_weighted, window, noverlap, nfft, fs);

% Plot PSDs
% subplot(3,1,1);
% semilogx(f, 10*log10(pxx_x), 'b', f, 10*log10(pxx_xw), 'r');
% title('Longitudinal Acceleration PSD');
% xlabel('Frequency (Hz)');
% ylabel('PSD (dB/Hz)');
% legend('Raw', 'Weighted');
% grid on;
% xlim([0.1 100]);
% 
% subplot(3,1,2);
% semilogx(f, 10*log10(pxx_y), 'b', f, 10*log10(pxx_yw), 'r');
% title('Lateral Acceleration PSD');
% xlabel('Frequency (Hz)');
% ylabel('PSD (dB/Hz)');
% legend('Raw', 'Weighted');
% grid on;
% xlim([0.1 100]);
% 
% subplot(3,1,3);
semilogx(f, 10*log10(pxx_z), 'b');
title('Vertical Acceleration PSD');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
legend('Raw', 'Weighted');
grid on;
xlim([0.1 100]);


hold on;

load flightData_AC.mat
load Dryden_opt_flightData_AC.mat
load flightData_AC_v0.mat
load flightData_UnOpt.mat
load flightData_AC.mat
% compute the ride quality of hte 
ax = flight_data.ax(:); % Ensure column vector
ay = flight_data.ay(:);
az = flight_data.az(:);
fs = flight_data.fs;
t = flight_data.t;


% Apply frequency weighting with proper filter syntax
% Define numerator & denominator for Wk
numWk = [1 32 1859 793.9 81.89];
denWk = [1 2261 7007 20640 79.96];
% Analog to discrete at sample rate fs
fs = flight_data.fs;  % e.g., 200 Hz
Ts = 1/fs;
[bdWk, adWk] = bilinear(numWk, denWk, fs);

% Define numerator & denominator for Wd
numWd = [0 14.55 6.026 7.714];     % leading zero for proper order
denWd = [1 51.65 47.66 15.02];
[bdWd, adWd] = bilinear(numWd, denWd, fs);

% Filter your acceleration signals
az_weighted = filter(bdWk, adWk, az);     % vertical axis
ax_weighted = filter(bdWd, adWd, ax);     % longitudinal
ay_weighted = filter(bdWd, adWd, ay);     % lateral

% Remove initial transient (1 second)
transient = round(1*fs);
t = t(transient:end);
ax_weighted = ax_weighted(transient:end);
ay_weighted = ay_weighted(transient:end);
az_weighted = az_weighted(transient:end);


% Calculate PSD using Welch's method
nfft = 2^nextpow2(length(t)/10);
window = hann(nfft);
noverlap = nfft/2;

[pxx_x, f] = pwelch(ax, window, noverlap, nfft, fs);
[pxx_xw, ~] = pwelch(ax_weighted, window, noverlap, nfft, fs);

[pxx_y, ~] = pwelch(ay, window, noverlap, nfft, fs);
[pxx_yw, ~] = pwelch(ay_weighted, window, noverlap, nfft, fs);

[pxx_z, ~] = pwelch(az, window, noverlap, nfft, fs);
[pxx_zw, ~] = pwelch(az_weighted, window, noverlap, nfft, fs);

% Plot PSDs
% subplot(3,1,1);
% semilogx(f, 10*log10(pxx_x), 'c', f, 10*log10(pxx_xw), 'm');
% title('Longitudinal Acceleration PSD');
% xlabel('Frequency (Hz)');
% ylabel('PSD (dB/Hz)');
% legend('Raw', 'Weighted');
% grid on;
% xlim([0.1 100]);
% 
% subplot(3,1,2);
% semilogx(f, 10*log10(pxx_y), 'c', f, 10*log10(pxx_yw), 'm');
% title('Lateral Acceleration PSD');
% xlabel('Frequency (Hz)');
% ylabel('PSD (dB/Hz)');
% legend('Raw', 'Weighted');
% grid on;
% xlim([0.1 100]);

% subplot(3,1,3);
semilogx(f, 10*log10(pxx_z), 'c');
title('Vertical Acceleration PSD');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
legend('Raw', 'Weighted');
grid on;
xlim([0.1 100]);

legend(["Rigid" "DAC"], "location","southwest")
