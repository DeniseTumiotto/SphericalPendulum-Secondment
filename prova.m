% main for frb
clearvars
close all
clc

%% physical parameters

% inertia tensor
% I1 = rand(1);
% I2 = rand(1);
% if I2 < I1
%     tmp = I1;
%     I1 = I2;
%     I2 = tmp;
% end
I1 = 0.2;
I2 = 0.5;
I3 = 1;

inertia = diag([1/I1, 1/I2, 1/I3]);

% initial angular momentum
m0 = 1 - rand(3, 1);
m0 = m0./norm(m0);

% attitude matrix
Q = eye(3);
Q = reshape(Q, 9, 1);

% damping value
alpha = 10;

%% numerical parameters

t0 = 0;
te = 50;
atol = 1e-10;
rtol = 1e-10;
% method = 'implicit Lie Euler method';
method = 'implicit midpoint rule';

%% initialization

% y0 = [m0; Q];
y0 = m0;

% LieSol = zeros(3+9, nTime);
% LieSol(:, 1) = y0;

%% equations and useful functions

% f = @(y) freeRigidBody(y, inertia);
f = @(y) -(eye(3)-y*transpose(y))*(inertia*y) - alpha*y;
fManiToAl = @(y) skw(y) * f(y);
action = @(a, b) a * b;
exponentialMap = @(x) expSO3(x);
res = @(yOld, yNew, h) residual(yOld, yNew, h, fManiToAl, action, exponentialMap, method);
jac = @(yOld, yNew, h) jacobian(yOld, yNew, h, fManiToAl, action, exponentialMap, method);

%% solution

options = odeset('AbsTol', atol, 'RelTol', rtol);
ySol = ode45(@(t, y) f(y), [t0, te], y0, options);

nTime = size(ySol.x, 2);
time = linspace(ySol.x(1), ySol.x(end), nTime);
zSol = ySol.y;
mSol = deval(ySol, time, 1:3);
% QSol = deval(ySol, time, 3 + (1:9));

LieSol = zeros(size(mSol));
LieSol(:, 1) = m0;
for i = 2:nTime
%     LieSol(:, i) = LieEuler(fManiToAl, action, exponentialMap, LieSol(:, i-1), time(i)-time(i-1));
    LieSol(:, i) = NewtonIt(res, jac, LieSol(:, i-1), time(i)-time(i-1), 100, atol, rtol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 2:nTime
%     LieSol(:, i) = LieEuler(f, action, exponentialMap, LieSol(:, i-1), dt);
% end

%% energy evaluation

kinEnergy = zeros(nTime, 1);
dEnergy = zeros(nTime, 1);
LieKinEnergy = zeros(nTime, 1);
LieDEnergy = zeros(nTime, 1);
energy = @(m) 0.5 * m' * inertia * m;
derenergy = @(m) -(inertia * m)' * (eye(3) - m*m')*(inertia*m); 

for i = 1:nTime
    kinEnergy(i) = energy(mSol(:, i));
    dEnergy(i) = derenergy(mSol(:, i));
    LieKinEnergy(i) = energy(LieSol(:, i));
    LieDEnergy(i) = derenergy(LieSol(:, i));
end

%% plots

figure()
plot(time, mSol(1, :), time, mSol(2, :), time, mSol(3, :), LineWidth=1.5)
legend('mx','my','mz')
grid on
title('Angular Momenta', 'FontSize', 20)
xlabel('time', FontSize=16)
ylabel('angular momenta', FontSize=16)

figure()
plot(time, LieSol(1, :), time, LieSol(2, :), time, LieSol(3, :), LineWidth=1.5)
legend('Lie mx', 'Lie my', 'Lie mz')
grid on
title('Angular Momenta', 'FontSize', 20)
xlabel('time', FontSize=16)
ylabel('angular momenta', FontSize=16)

figure()
plot(time, vecnorm(mSol), time, vecnorm(LieSol), LineWidth=1.5)
legend('sol norm', 'Lie sol norm')
grid on
title('Norm of the solution', 'FontSize', 20)
xlabel('time', FontSize=16)
ylabel('angular momenta', FontSize=16)

figure()
plot(time, kinEnergy, LineWidth=1.5)
hold on
plot(time, LieKinEnergy, LineWidth=1.5)
grid on
legend('Energy','Lie Energy')
title('Kinetic Energy', 'FontSize', 20)
xlabel('time', FontSize=16)
ylabel('Kinetic energy', FontSize=16)

figure()
plot(time, dEnergy, LineWidth=1.5)
hold on
plot(time, LieDEnergy, LineWidth=1.5)
grid on
legend('Energy derivative','Lie Energy derivative')
title('Kinetic Energy', 'FontSize', 20)
xlabel('time', FontSize=16)
ylabel('Kinetic energy', FontSize=16)
