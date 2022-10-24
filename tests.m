%% number of tests

N = 2;

%% METHOD
% method = 1 --> explicit Lie Euler
% method = 2 --> implicit Lie Euler
% method = 3 --> implicit Lie midpoint rule
% method = 4 --> ode45 (MATLAB routine)
% method = 5 --> closes run
% method > 5 --> gives error

% method = [1 1 1 1 1 1];
method = 1*ones(N,1);

%% TIME and STEPS
% T = end time of the simulation
% N_TIME = number of steps
T = 10;
N_TIME = 1e6 * ones(1, N);
% T = 1;
% N_TIME = logspace(2, 5, N-1);
% N_TIME = ceil(N_TIME);
% N_TIME(end+1) = 1e6;

%% system parameters

k = 0;
% damp = linspace(10, 1000, N-5);
damp = 1 * ones(N, 1);
stiff = 1e2 * ones(N, 1);

%% initial conditions

% [~, ~, z0(:, 1)] = initializeZeroVel();
% [~, ~, z0(:, 2)] = initializeSmallVariation(z0(:, 1));
z0 = zeros(6, N);
% [~, ~, z0(:, 1)] = initializeZeroVel();
z0(:, 1) = startingValue();
z0(:, 2) = startingValue();
% for i = 2:N
%     z0(:, i) = z0(:, 1);
% end
% [~, ~, z0(:, 2)] = initializeSmallVariation(z0(:, 1), 0, 0.1);

for i = 1:N-2
    z0(:, 2*(i-1)+1) = z0(:, 1);
    z0(:, 2*i) = z0(:, 2);
end

for i = 1:N
    disp(['(', num2str(i), ')'])
    main(method(i), N_TIME(i), k, damp(i), stiff(i), z0(:, i), T)
    pause(1)
end


clearvars
