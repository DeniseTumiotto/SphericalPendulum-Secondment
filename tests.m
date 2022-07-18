%% number of tests

N = 6;

%% METHOD
% method = 1 --> explicit Lie Euler
% method = 2 --> implicit Lie Euler
% method = 3 --> implicit Lie midpoint rule
% method = 4 --> ode45 (MATLAB routine)
% method = 5 --> closes run
% method > 5 --> gives error

method = [1 1 2 2 1 1];

%% TIME and STEPS
% T = end time of the simulation
% N_TIME = number of steps
T = 10;
N_TIME = 1e5 * ones(1, N);
% T = 1;
% N_TIME = logspace(0.01, 3, N);
% N_TIME = ceil(N_TIME);
% N_TIME(end+1) = 1e6;

%% system parameters

k = 0;
% damp = linspace(10, 1000, N-5);
damp = [10 10 10 10 1 1];

%% initial conditions

% [~, ~, z0(:, 1)] = initializeZeroVel();
% [~, ~, z0(:, 2)] = initializeSmallVariation(z0(:, 1));
z0 = zeros(6, N);
[~, ~, z0(:, 1)] = initializeZeroVel();
[~, ~, z0(:, 2)] = initializeSmallVariation(z0(:, 1), 0, 0.1);

for i = 1:4
    z0(:, 2*(i-1)+1) = z0(:, 1);
    z0(:, 2*i) = z0(:, 2);
end

for i = 1:N
    main(method(i), N_TIME(i), k, damp(i), z0(:, i), T)
    pause(1)
end


clearvars
