% parameters

lambda = -0.01;
mu = -1;

eq = @(t, y) [lambda * y(1); mu * y(2)];

t0 = 0;
te = 1;
nStep = 1000;
timeVec = linspace(t0, te, nStep);

% initial condition

y0 = [0; -pi/2];

% solution

[tSol, ySol] = ode45(eq, timeVec, y0);

% plot

figure('Units','centimeters','Position',[10 10 20 15])
[xS2, yS2, zS2] = sphere(360);
h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
h.EdgeColor = 'none';
hold on
for i = 1:nStep
    % to cart coord
    q = sph2cart(ySol(i, :));
    % Plot pendulum at time step i
    plot3(q(1), q(2), q(3), '*r', 'MarkerSize', 2);
    xlabel("x")
    ylabel("y")
    zlabel("z")
    str = "Time evolution, T="+num2str((timeVec(2)-timeVec(1))*i);
    axis equal;
    title(str, 'FontSize', 16)
    pause(0.0000000000001);
end
