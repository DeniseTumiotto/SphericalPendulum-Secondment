syms x y z
syms phi theta
syms lam mu


theta = atan(y/x);
phi = asin(z);

vecField = cross([sin(theta)*lam*phi; -cos(theta)*lam*phi; mu*theta], [x; y; z]);

assume(mu<0)
assume(lam<0)

assume(x == -sqrt(2)/2)
assume(y == 0)
assume(z == sqrt(2)/2)

disp(sign(vecField))

figure()
% Create sphere surface
[xS2, yS2, zS2] = sphere(360);
h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
h.EdgeColor = 'none';
hold on
plot3(sqrt(2)/2,0,sqrt(2)/2, 'or', 'MarkerSize',5, 'MarkerFaceColor', 'r')
plot3(1,0,0, 'ok', 'MarkerSize',5, 'MarkerFaceColor', 'k')
xlabel('x')
ylabel('y')
zlabel('z')