% main for frb
clearvars
close all

%% physical parameters

% inertia tensor
I1 = rand(1);
I2 = rand(1);
if I2 > I1
    tmp = I1;
    I1 = I2;
    I2 = tmp;
end
I3 = 1;

inertia = diag([1/I1, 1/I2, 1/I3]);

% initial angular momentum
m0 = rand(3, 1);
m0 = m0./norm(m0);

% attitude matrix
Q = eye(3);
Q = reshape(Q, 9, 1);

% damping value
% alpha = 0.1;

%% numerical parameters

t0 = 0;
te = 100;
atol = 1e-13;
rtol = 1e-13;

%% initialization

y0 = [m0; Q];

% LieSol = zeros(3+9, nTime);
% LieSol(:, 1) = y0;

%% equations and useful functions

if exist('alpha','var')
    f = @(y) freeRigidBody(y, inertia, alpha);
else
    f = @(y) freeRigidBody(y, inertia);
end
action = @(a, b) actionSO3xR3(a, b);
exponentialMap = @(x) expSO3xR3(x);

%% solution

options = odeset('AbsTol', atol, 'RelTol', rtol);
ySol = ode45(@(t, y) f(y), [t0, te], y0);

nTime = size(ySol.x, 2);
time = linspace(ySol.x(1), ySol.x(end), nTime);
zSol = ySol.y;
mSol = deval(ySol, time, 1:3);
QSol = deval(ySol, time, 3 + (1:9));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 2:nTime
%     LieSol(:, i) = LieEuler(f, action, exponentialMap, LieSol(:, i-1), dt);
% end

%% energy evaluation

kinEnergy = zeros(nTime, 1);
energy = @(m) 0.5 * m' * inertia * m;

for i = 1:nTime
    kinEnergy(i) = energy(mSol(:, i));
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
plot(time, kinEnergy, LineWidth=1.5)
grid on
title('Kinetic Energy', 'FontSize', 20)
xlabel('time', FontSize=16)
ylabel('Kinetic energy', FontSize=16)


% Qm = zeros(3, nTime);
% for i = 1:nTime
%     Qm(:, i) = reshape(QSol(:, 1), 3, 3) * mSol(:, i);
% end
% figure()
% plot(time, Qm(1, :), time, Qm(2, :), time, Qm(3, :), LineWidth=1.5)
% title('Spatial Angular Momentum', "FontSize", 20)
% grid on
% xlabel('time', FontSize=16)
% ylabel('spatial angular momenta', FontSize=16)

% figure()
% plot(time, vecnorm(mSol), LineWidth=1.5)
% grid on
% title('Norm of Angular Momentum', 'FontSize', 20)
% xlabel('time', FontSize=16)
% ylabel('angular momenta', FontSize=16)

%% animation

% animation(QSol, time)
% 
% function animation(sols, time)
% e1 = [1; 0; 0];
% 
% figure('Units','centimeters','Position',[10 10 20 15])
% % Create sphere surface
% [xS2, yS2, zS2] = sphere(360);
% h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
% h.EdgeColor = 'none';
% hold on
% crr_rot = reshape(sols(:, 1), 3, 3) * e1;
% plot3(crr_rot(1), crr_rot(2), crr_rot(3), 'o');
% for i = 2:size(time, 2)
%     old_rot = crr_rot;
%     % Plot pendulum at time step i
%     crr_rot = reshape(sols(:, i), 3, 3) * e1;
%     plot3([old_rot(1) crr_rot(1)], [old_rot(2) crr_rot(2)], [old_rot(3) crr_rot(3)], '-');
%     xlabel("x")
%     ylabel("y")
%     zlabel("z")
%     str = "Time evolution of the rotation";
%     axis equal;
%     title(str, 'FontSize', 16)
%     grid on
%     pause(0.000000001);
% end
% end
