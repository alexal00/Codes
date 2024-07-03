% Generate Lorenz attractor
sigma = 16;
rho = 28;
beta = 8;
x0 = [2; 3; -14];
tspan = [0 100];
[t, sol] = solveLorenzODE(sigma, rho, beta, x0, tspan);

% Extract DOFs
[x, y, z] = extractDOFs(sol);

% Define parameters for MLCE calculation
max_dim = 15;
max_tau = 100;
Rtol = 10;
Atol = 2;
fs = 1 / (t(2) - t(1));  % Sampling frequency

% Call calculateMLCE function
calculateMLCE([x y z], max_dim, max_tau, Rtol, Atol, fs);

% Function to solve Lorenz system of ODEs
function [t, sol] = solveLorenzODE(sigma, rho, beta, x0, tspan)
    lorenz = @(t, x) [sigma * (x(2) - x(1));
                      x(1) * (rho - x(3)) - x(2);
                      x(1) * x(2) - beta * x(3)];
    [t, sol] = ode45(lorenz, tspan, x0);
end

% Function to extract DOFs
function [x, y, z] = extractDOFs(sol)
    x = sol(:, 1);
    y = sol(:, 2);
    z = sol(:, 3);
end


function calculateMLCE(timeSeries, max_dim, max_tau, Rtol, Atol, fs)
    m = estimateEmbeddingDimension(timeSeries, max_dim, max_tau, Rtol, Atol);
    tau = estimateDelay(timeSeries, max_tau);
    X = reconstructPhaseSpace(timeSeries, m, tau);
    Tmean = meanPeriod(X, fs);
    neighbors = findNearestNeighbors(X, Tmean);
    lyapunovExponent = estimateLyapunovExponent(X, neighbors, tau, fs);
    fprintf('Optimal Embedding Dimension: %d\n', m);
    fprintf('Optimal Time Delay (tau): %d\n', tau);
    fprintf('Estimated Maximum Lyapunov Exponent: %f\n', lyapunovExponent);
end

function m = estimateEmbeddingDimension(timeSeries, max_dim, max_tau, Rtol, Atol)
    fnnPercentages = zeros(1, max_dim);
    for d = 1:max_dim
        fnnPercentages(d) = fnn(timeSeries, d, max_tau, Rtol, Atol);
    end
    m = find(fnnPercentages < 1, 1);
    if isempty(m)
        m = max_dim;
    end
end

function fnnPercentage = fnn(timeSeries, m, max_tau, Rtol, Atol)
    timeSeries = reshape(timeSeries, [], size(timeSeries, 3)); % Reshape timeSeries if needed
    N = size(timeSeries, 1) - (m+1) * max_tau;
    Y = zeros(N, m);
    for i = 1:m
        Y(:, i) = timeSeries((1:N) + (i-1) * max_tau);
    end

    Z = timeSeries((1:N) + m * max_tau);

    fnnCount = 0;
    totalCount = 0;

    for i = 1:N
        distances = vecnorm(Y - Y(i, :), 2, 2);
        [sortedDistances, indices] = sort(distances);

        r1 = sortedDistances(2);
        nearestNeighbor = indices(2);

        r2 = sqrt(r1^2 + (Z(i) - Z(nearestNeighbor))^2);

        if (r2 / r1 > Rtol) || (abs(Z(i) - Z(nearestNeighbor)) / max(std(timeSeries(:)), eps) > Atol)
            fnnCount = fnnCount + 1;
        end

        totalCount = totalCount + 1;
    end

    fnnPercentage = (fnnCount / totalCount) * 100;
end

function tau = estimateDelay(timeSeries, max_tau)
    amiValues = zeros(1, max_tau);
    for lag = 1:max_tau
        amiValues(lag) = average_mutual_information(timeSeries,lag);
    end
    [~, tau] = min(amiValues);
end

function ami = average_mutual_information(x, lag)
    % Compute average mutual information for a given lag
    x1 = x(1:end-lag);
    x2 = x(lag+1:end);

    pX1 = histcounts(x1, 'Normalization', 'probability');
    pX2 = histcounts(x2, 'Normalization', 'probability');
    pX1X2 = histcounts2(x1, x2, 'Normalization', 'probability');

    pX1 = pX1(pX1 > 0);
    pX2 = pX2(pX2 > 0);
    pX1X2 = pX1X2(pX1X2 > 0);

    Hx1 = -sum(pX1 .* log2(pX1));
    Hx2 = -sum(pX2 .* log2(pX2));
    Hx1x2 = -sum(pX1X2 .* log2(pX1X2));

    ami = Hx1 + Hx2 - Hx1x2;
end

function X = reconstructPhaseSpace(timeSeries, m, tau)
    N = size(timeSeries, 1) - (m-1) * tau;
    X = zeros(N, m);
    for i = 1:m
        X(:, i) = timeSeries((1:N) + (i-1) * tau);
    end
end

function Tmean = meanPeriod(X, fs)
    [pxx, f] = periodogram(X(:, 1), [], [], fs);
    [~, maxIdx] = max(pxx);
    dominantFreq = f(maxIdx);
    Tmean = ceil(fs / dominantFreq);
end

function neighbors = findNearestNeighbors(X, Tmean)
    N = size(X, 1);
    neighbors = zeros(N, 1);
    for i = 1:N
        distances = vecnorm(X - X(i, :), 2, 2);
        validIndices = find(abs((1:N) - i) > Tmean);
        if ~isempty(validIndices)
            [~, minIdx] = min(distances(validIndices));
            neighbors(i) = validIndices(minIdx);
        else
            neighbors(i) = NaN;
        end
    end
end

function lyapunovExponent = estimateLyapunovExponent(X, neighbors, tau, fs)
    N = size(X, 1);
    validIndices = ~isnan(neighbors);
    divergence = zeros(sum(validIndices), N - tau);
    validNeighbors = neighbors(validIndices);
    for i = 1:sum(validIndices)
        for k = 1:N - tau
            if validNeighbors(i) + k <= N && i + k <= N
                divergence(i, k) = log(norm(X(i + k, :) - X(validNeighbors(i) + k, :)));
            end
        end
    end
    lyapunovExponent = mean(divergence(:)) / tau;
end
