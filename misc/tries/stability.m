% Load or generate your time history data
% Assume t is time vector and x is the time history data
t = linspace(0, 100, 1000); % Example time vector
x = 2*exp(-0.1*t).*sin(t) ;%+ 0.1 * randn(size(t)); % Example data (sinusoidal + noise)

% Plot the time history
figure;
plot(t, x);
title('Time History');
xlabel('Time');
ylabel('Amplitude');

% Statistical analysis
mean_x = mean(x);
std_x = std(x);
fprintf('Mean: %f\n', mean_x);
fprintf('Standard Deviation: %f\n', std_x);

% Frequency analysis using FFT
n = length(x);
f = (0:n-1)*(1/(t(2)-t(1)))/n; % Frequency vector
X = fft(x);
amplitude = abs(X)/n;

figure;
plot(f, amplitude);
title('Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Phase space plot (x vs. dx/dt)
dx_dt = gradient(x, t);
figure;
plot(x, dx_dt);
title('Phase Space');
xlabel('x');
ylabel('dx/dt');

% Lyapunov Exponent (Simplified calculation)
maxIterations = length(x);
d0 = 1e-4; % Initial separation
lyapunovExponent = 0;
for i = 1:maxIterations-1
    d1 = abs(x(i+1) - x(i));
    if d1 > d0
        lyapunovExponent = lyapunovExponent + log(d1/d0);
    end
end
lyapunovExponent = lyapunovExponent / (maxIterations - 1);

fprintf('Lyapunov Exponent: %f\n', lyapunovExponent);

% Interpretation
if lyapunovExponent < 0
    disp('The system is stable.');
elseif lyapunovExponent == 0
    disp('The system may be leading to an LCO.');
else
    disp('The system is unstable.');
end
