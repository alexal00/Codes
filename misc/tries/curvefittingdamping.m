close all; clear variables; clc
% Generate example data: simulated multi-harmonic damped sinusoidal signal
time = linspace(0, 10, 1000);
A0 = 1.0;  % Initial amplitude
zeta_true = 0.05;  % True damping ratio
omega_n = 2 * pi;  % Fundamental natural frequency (rad/s)
harmonics = [1, 2, 3];  % Harmonics as integer multiples of the fundamental frequency

signal = zeros(size(time));
for k = harmonics
    signal = signal + A0 * exp(-zeta_true * omega_n * time) .* cos(k * omega_n * sqrt(1 - zeta_true^2) * time);
end

% Add some noise to the signal
noise = 0.05 * randn(size(time));
signal_noisy = signal + noise;

% Plot the noisy signal
figure;
plot(time, signal_noisy);
xlabel('Time (s)');
ylabel('Amplitude');
title('Noisy Multi-Harmonic Damped Signal');
grid on;

% Hilbert Transform Method
analytic_signal = hilbert(signal_noisy);
envelope = abs(analytic_signal);

% Fit an exponential decay function to the envelope
fit_func = @(b, t) b(1) * exp(-b(2) * t);  % Exponential decay function
initial_guess = [max(envelope), 0.1];  % Initial guess for parameters

% Define the objective function for fitting
objective_func = @(b) sum((envelope - fit_func(b, time)).^2);

% Perform non-linear least squares fitting
opts = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8, 'MaxFunEvals', 1000, 'MaxIter', 1000);
beta = fminsearch(objective_func, initial_guess, opts);

% Extract estimated damping ratio and amplitude
A0_fit = beta(1);
lambda_fit = beta(2);

% Damping ratio estimation
zeta_estimated_hilbert = lambda_fit / omega_n;

disp(['Hilbert Method - Estimated Initial Amplitude: ', num2str(A0_fit)]);
disp(['Hilbert Method - Estimated Damping Ratio: ', num2str(zeta_estimated_hilbert)]);

% Plot fitted envelope
hold on;
plot(time, fit_func(beta, time), 'g--', 'LineWidth', 2);
legend('Noisy Multi-Harmonic Damped Signal', 'Envelope', 'Fitted Envelope');
hold off;

% Logarithmic Decrement Method
[peaks, locs] = findpeaks(signal_noisy, time,'MinPeakDistance',1);

% Calculate logarithmic decrement
num_peaks = length(peaks);
delta = 1 / (num_peaks - 1) * sum(log(peaks(1:end-1) ./ peaks(2:end)));

% Damping ratio estimation
zeta_estimated_logdec = delta / sqrt(4 * pi^2 + delta^2);

disp(['Logarithmic Decrement Method - Estimated Damping Ratio: ', num2str(zeta_estimated_logdec)]);



