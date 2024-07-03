% Generate Lorenz attractor data
dt = 0.01;
t = 0:dt:50;
sigma = 10;
beta = 8/3;
rho = 28;
x0 = [0; 1; 20];
f = @(t,x) [sigma*(x(2)-x(1)); x(1)*(rho-x(3))-x(2); x(1)*x(2)-beta*x(3)];
[~,x] = ode45(f, t, x0);

% Save Lorenz attractor data
save('lorenz_attractor_data.mat', 't', 'x');

% Calculate maximum Lyapunov exponent using Rosenstein method
max_lyap_exp = rosenstein_method(x(:,1), 1/dt);

% Display the result
disp(['Estimated Maximum Lyapunov Characteristic Exponent: ', num2str(max_lyap_exp)]);

% False Nearest Neighbors algorithm
function [embedding_dim, tau] = f_nn(x, m_max, tau_max)
    % Initialize parameters
    N = length(x);
    n_max = min(2000, N-1); % Maximum number of points to consider
    
    % Initialize arrays to store results
    epsilon = zeros(1, tau_max);
    R = zeros(1, tau_max);
    n_nn = zeros(1, tau_max);
    
    % Calculate phase space distance
    for tau = 1:tau_max
        phase_space = zeros(N-(m_max-1)*tau, m_max);
        for i = 1:m_max
            phase_space(:,i) = x((i-1)*tau+1:N-(m_max-i)*tau);
        end
        
        for i = 1:N-(m_max-1)*tau
            point = phase_space(i,:);
            distances = sqrt(sum((phase_space - point).^2,2));
            [~,idx] = sort(distances);
            nn = idx(2);
            % Check if indices are within bounds
            if i+m_max*tau <= N && nn+m_max*tau <= N
                if abs(x(i+m_max*tau) - x(nn+m_max*tau)) > 0
                    epsilon(tau) = epsilon(tau) + distances(nn);
                    R(tau) = R(tau) + abs(x(i+m_max*tau) - x(nn+m_max*tau));
                    n_nn(tau) = n_nn(tau) + 1;
                end
            end
        end
        epsilon(tau) = epsilon(tau) / n_nn(tau);
        R(tau) = R(tau) / n_nn(tau);
    end
    
    % Find the first minimum of R/epsilon
    R_over_epsilon = R ./ epsilon;
    [~, tau] = min(R_over_epsilon);
    
    % Find the minimum embedding dimension using mutual information
    mi = zeros(1, m_max);
    for m = 1:m_max
        phase_space = zeros(N-(m-1)*tau, m);
        for i = 1:m
            phase_space(:,i) = x((i-1)*tau+1:N-(m-i)*tau);
        end
        mi(m) = mutual_information(phase_space);
    end
    embedding_dim = find(mi(1:end-1) > mi(2:end), 1) + 1;
end


% Mutual information estimation
function mi = mutual_information(X)
    [n, m] = size(X);
    edges = cell(1, m);
    for i = 1:m
        edges{i} = linspace(min(X(:,i)), max(X(:,i)), 10);
    end
    pxy = hist3(X, 'Edges', edges) / n;
    px = sum(pxy, 2);
    py = sum(pxy, 1);
    px_py = px * py;
    pxy(pxy == 0) = 1; % Avoid log(0) issues
    mi = sum(sum(pxy .* log2(pxy ./ px_py)));
end

% Rosenstein method for calculating maximum Lyapunov exponent
function max_lyap_exp = rosenstein_method(x, fs)
    % Parameters
    m_max = 10; % Maximum embedding dimension to search for
    tau_max = 100; % Maximum delay to search for

    % Calculate False Nearest Neighbors
    [embedding_dim, tau] = f_nn(x, m_max, tau_max);

    % Construct trajectory matrix
    X = trajectory_matrix(x, embedding_dim, tau);

    % Calculate nearest neighbor threshold
    mean_period = ceil(fs / max(diff(t)));
    nearest_neighbor_threshold = mean_period * 2; % Adjust the multiplier as needed

    % Calculate nearest neighbors
    [m, N] = size(X);
    distances = zeros(N, N);
    for i = 1:N
        for j = 1:N
            distances(i, j) = norm(X(:, i) - X(:, j));
        end
    end

    % Find nearest neighbors with temporal separation greater than mean_period
    nearest_neighbors = zeros(N, 1);
    for i = 1:N
        temp = distances(i, :);
        temp(temp == 0) = NaN; % Exclude self-distances
        [sorted_distances, ~] = sort(temp);
        idx = find(sorted_distances > nearest_neighbor_threshold, 1);
        nearest_neighbors(i) = idx;
    end

    % Calculate the maximum Lyapunov exponent
    max_lyap_exp = lyapunov_exponent(X, nearest_neighbors);
end

% Trajectory matrix construction
function X = trajectory_matrix(x, m, tau)
    N = length(x);
    X = zeros(m, N - (m - 1) * tau);
    for i = 1:m
        X(i, :) = x((i - 1) * tau + 1:N - (m - i) * tau);
    end
end

% Lyapunov exponent calculation
function max_lyap_exp = lyapunov_exponent(X, nearest_neighbors)
    [m, N] = size(X);
    max_lyap_exp = 0;
    for i = 1:N
        for j = i + 1:N
            dist = norm(X(:, i) - X(:, j));
            if dist > nearest_neighbors(i)
                delta = dist / norm(X(:, i) - X(:, j));
                max_lyap_exp = max(max_lyap_exp, log(delta));
            end
        end
    end
    max_lyap_exp = max_lyap_exp / (N * (N - 1) / 2);
end
