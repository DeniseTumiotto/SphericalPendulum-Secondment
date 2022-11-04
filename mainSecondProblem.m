h = 0.1:0.1:15;
[~, N] = size(h);

tol = 1e-10;

mu = -2;
lam = -0.1;
% phi0 = pi/2*(2*rand(2,1)-1);
% theta0 = 2*pi*rand(2,1);
phi0 = [0, 0];
theta0 = [pi/2+0.1, pi/2-0.1];
y01 = sph2cart([phi0(1);theta0(1)]);
y02 = sph2cart([phi0(2);theta0(2)]);

y1 = zeros(3,N+1);
y2 = zeros(3,N+1);
y1(:,1) = y01;
y2(:,1) = y02;

yi1 = zeros(3,N+1);
yi2 = zeros(3,N+1);
yi1(:,1) = y01;
yi2(:,1) = y02;

A = @(y) matrix(y,lam,mu);

lieStep = @(y0,y,dt) expm(dt*A(y))*y0;

implicitLie = @(y0,y,dt) -y+lieStep(y0,y,dt);

for index = 1:N
    y1(:,index+1) = lieStep(y1(:,1),y1(:,1),h(index));
    y2(:,index+1) = lieStep(y2(:,1),y2(:,1),h(index));
    
    yi1(:,index+1) = fsolve(@(x) implicitLie(yi1(:,1),x,h(index)),yi1(:,index));
    yi2(:,index+1) = fsolve(@(x) implicitLie(yi2(:,1),x,h(index)),yi2(:,index));
end

plotSecondProblem(yi1,yi2)

function rslt = matrix(y,lam,mu)
    s = cart2sph(y);
    rslt = skw([sin(s(3))*lam*s(2); -cos(s(3))*lam*s(2); mu*s(3)]);
end