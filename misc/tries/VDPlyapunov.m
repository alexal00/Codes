% Define parameters for the Van der Pol oscillator
mu = 1;  % parameter controlling the nonlinearity of the oscillator

% Define the Van der Pol system of ODEs
van_der_pol = @(t, y) [y(2); mu*(1 - y(1)^2)*y(2) - y(1)];

% Set initial conditions
y0 = [0.1; 0.1];  % initial state

% Time span for simulation
tspan = [0 100];

% Solve the system numerically
[t, y] = ode45(van_der_pol, tspan, y0);

% Define a small perturbation for each component of the state
epsilon = 1e-6;
perturbation = epsilon * eye(2);

% Initialize the matrix to store Lyapunov exponents
num_steps = length(t);
num_exponents = size(perturbation, 1);
exponents = zeros(num_exponents, 1);

% Iterate over each time step to compute the Lyapunov exponents
for i = 1:num_steps
    % Jacobian matrix of the system at current state
    J = [0, 1; -2*mu*y(i, 1)*y(i, 2) - 1, mu*(1 - y(i, 1)^2)];

    % Perform QR decomposition to orthogonalize the perturbation matrix
    [~, R] = qr(J * perturbation);
    
    % Update the perturbation matrix
    perturbation = R;
    
    % Renormalize the perturbation matrix
    perturbation = perturbation / norm(perturbation(1, :));
    
    % Compute the accumulated logarithm of the diagonal elements of perturbation matrix
    exponents = exponents + log(abs(diag(R)));
end

% Average over time to get the Lyapunov exponents
exponents = exponents / num_steps;

% Display the Lyapunov exponents
disp('Lyapunov Exponents:');
disp(exponents);
figure
plot(y(:,1),y(:,2))