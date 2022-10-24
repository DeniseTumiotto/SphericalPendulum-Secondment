clearvars
close all
clc

% physical parameters
d = 1;
k = 1;
m = 1;

% time intergration parameters
t0 = 0;
te = 10;
nTime = 1e5;
dt = (te-t0)/nTime;
timeVec = linspace(t0, te, nTime);

% save parameters
timestamp = string(datetime('now','format','yyyyMMdd''T''HHmmss'));
filename = strcat('out/', timestamp, 'lin', 'prm', '.mat');
save(filename)

% initial condition
[~, ~, y0] = initializeZeroVel();
txt = input('small variation instead? [Y/N]', "s");
if strcmp(txt, "Y") || strcmp(txt, "")
    [~, ~, y0] = initializeSmallVariation(startingValue(), false, 0.1);
end


% equations of motion
eq = @(y) [     zeros(3)         eye(3);
           -k/m*eye(3)      -d/m*eye(3)] * y;

% equations of motion
ode = @(t, y) [     zeros(3)         eye(3);
               -k/m*eye(3)      -d/m*eye(3)] * y;

action = @(B, input) actionSE3(B, input);

% solve with forward Euler
% zSol = zeros(6, nTime);
% zSol(:, 1) = y0;
% for i = 2:nTime
%     zSol(:, i) = LieEulerSE3(eq, action, zSol(:, i-1), dt);
% end

% solve with ode45
[tout, zSol] = ode45(ode, timeVec, y0);
zSol = transpose(zSol);

% save solution
filename = strcat('out/', timestamp, 'lin', 'sol', '.mat');
save(filename, 'zSol')
