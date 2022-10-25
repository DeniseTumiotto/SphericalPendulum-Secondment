% parameters

lambda = -0.001;
mu = -10;

eq = @(t, y) [lambda * y(1); mu * y(2)];

t0 = 0;
te = 1000;
nStep = 100000;
timeVec = linspace(t0, te, nStep);
dt = timeVec(2)-timeVec(1);

% save parameters

timestamp = string(datetime('now','format','yyyyMMdd''T''HHmmss'));
filename = strcat('out/', timestamp, 'S2', 'prm', '.mat');
save(filename)

% initial condition

y0 = 2*pi .* rand(2, 1) - pi;
% y0 = [1.5; -0.6];
zSol = zeros(2, nStep);
zSol(:, 1) = y0;

% functions

vecField = @(x) [lambda 0; 0 mu];
expMap = @(x) expm(x);
action = @(a, b) a*b;

% solution

% [tSol, ySol] = ode45(eq, timeVec, y0);
for i = 1:nStep
    zSol(:, i+1) = LieEuler(vecField, action, expMap, zSol(:, i), dt);
end

% save solution
filename = strcat('out/', timestamp, 'S2', 'sol', '.mat');
save(filename, 'zSol')

% plot

% figure('Units','centimeters','Position',[10 10 20 15])
% [xS2, yS2, zS2] = sphere(360);
% h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
% h.EdgeColor = 'none';
% hold on
% plot3(1, 0, 0, 'ob', 'MarkerSize', 3, 'MarkerEdgeColor','none', 'MarkerFaceColor', 'b')
% for i = 1:nStep/100:nStep
%     % to cart coord
%     q = sph2cart(zSol(:, i));
%     % Plot pendulum at time step i
%     plot3(q(1), q(2), q(3), '*r', 'MarkerSize', 2);
%     xlabel("x")
%     ylabel("y")
%     zlabel("z")
%     str = "Time evolution, T="+num2str(dt*i);
%     axis equal;
%     title(str, 'FontSize', 16)
%     pause(0.00000000000000001)
% end
