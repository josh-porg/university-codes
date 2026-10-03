% Generate random noise
n_samples = 1000;
noise = randn(1, n_samples);

% Apply a moving average filter to smooth the noise
window_size = 50; % Size of the moving window
smooth_noise = filter(ones(1, window_size) / window_size, 1, noise);
% filiter again on the fine scale to make it more differentiable
window_size = 5; % Size of the moving window
smooth_noise = filter(ones(1, window_size) / window_size, 1, smooth_noise);

% Generate time vector
time = linspace(0, 10, n_samples);

% Plot the original and smoothed noise
figure;
subplot(2, 1, 1);
plot(time, noise);
title('Original Random Noise');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(time, smooth_noise);
title('Smoothed Random Noise');
xlabel('Time');
ylabel('Amplitude');
