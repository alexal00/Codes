close all; clear all; clc
% Example data: simulated damped sinusoidal signal
time = linspace(0, 10, 1000);
A0 = 1.0;  % Initial amplitude
zeta_true = 0.05;  % True damping ratio
omega_n = 2 * pi;  % Natural frequency (rad/s)
signal = A0 * exp(-zeta_true * omega_n * time) .* cos(omega_n * sqrt(1 - zeta_true^2) * time);

% Add some noise to the signal
noise = 0.05 * randn(size(time));
signal_noisy = signal + noise;

% Plot the noisy signal
figure;
plot(time, signal_noisy);
xlabel('Time (s)');
ylabel('Amplitude');
title('Noisy Damped Signal');
grid on;

% Perform Fourier Transform
Fs = 1 / (time(2) - time(1));  % Sampling frequency
n = length(time);
Y = fft(signal_noisy);
f = (0:n-1) * (Fs / n);  % Frequency vector
amplitude_spectrum = abs(Y) / n;

% Plot amplitude spectrum
figure;
plot(f, amplitude_spectrum);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum');
grid on;

% Identify resonant peak and bandwidth
[~, peak_index] = max(amplitude_spectrum);
omega_r = f(peak_index) * 2 * pi;  % Resonant frequency in rad/s

% Find the -3dB points (half-power points)
half_power_amplitude = max(amplitude_spectrum) / sqrt(2);
above_half_power = amplitude_spectrum >= half_power_amplitude;

% Find frequency indices corresponding to half-power points
indices = find(above_half_power);
bandwidth = (f(indices(end)) - f(indices(1))) * 2 * pi;  % Bandwidth in rad/s

% Calculate damping ratio
zeta_estimated = bandwidth / (2 * omega_r);

disp(['Estimated Damping Ratio: ', num2str(zeta_estimated)]);
