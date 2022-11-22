syms phi theta t
syms a c

x = cos(phi(t))*cos(theta(t));
y = cos(phi(t))*sin(theta(t));
z = sin(phi);

m = @(xx) (a*cos(phi)^2+c*sin(phi)^2-xx);

dx = m(a)*cos(phi)*cos(theta);
dy = m(a)*cos(phi)*sin(theta);
dz = m(c)*sin(phi);
