h = 0.1:0.1:15;
[~, N] = size(h);

tol = 1e-10;

mu = -2;
lam = -0.1;
phi0 = pi/2*(2*rand(2,1)-1);
theta0 = 2*pi*rand(2,1);

A = @(y) matrix(y);

lieStep = @(y0,y,dt) expm(dt*A(y))*y0;

implicitLie = @(y0,y,dt) -y+lieStep(y0,y,dt);






function rslt = matrix(y)
s = cart2sph(y);
rslt = skw([sin(s(3))*lam*s(2); -cos(s(3))*lam*s(2); mu*s(3)]);
end