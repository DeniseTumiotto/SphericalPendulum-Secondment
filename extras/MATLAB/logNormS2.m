clearvars
clc

syms phi theta t
syms a c

x = cos(phi)*cos(theta);
y = cos(phi)*sin(theta);
z = sin(phi);

dphi = (a-c)*z*sqrt(x^2+y^2);
dtheta = (2*x*y)/(x^2+y^2) * (a*x+a*y+c*z-a);


B = [diff(dphi,phi) diff(dphi,theta);
     diff(dtheta,phi) diff(dtheta,theta)];

A = [0 dtheta*sin(phi)*cos(theta);
     -dtheta*tan(phi) -dphi*tan(phi)];

L = A+B;

G = [1 0;
     0 cos(phi)^2];

J = 0.5*((G*L)/G+transpose(L));

[V,D] = eig(J);

D = subs(D, {phi, theta}, {atan(y/x),asin(z)});
disp(simplify(diag(D)))
