syms a c
syms qx qy qz

h = a*qx^2 + a*qy^2 + c*qz^2;
hx = 2*a*qx;
hy = 2*a*qy;
hz = 2*c*qz;

A = [h+qx*hx-a 2*a*qx*qy (a+c)*qx*qz;
     2*qx*qy   h+qy*hy-a (a+c)*qy*qz;
     (a+c)*qx*qz (a+c)*qy*qz h+qz*hz-c];

lam = eig(A);

% assume(qx^2+qy^2+qz^2==1)
% assume(c<a)
disp(lam)