function F = fManiToAlgebra(y, damp, k, stiff)

% Assembling RHS of the mathematical pendulum
% subject to Earth gravitation [g = 9.81 m/s^2]
% Introducing some damping in the equations
q = y(1:3);
w = y(3 + (1:3));

g = 9.81;
e3 = [0; 0; 1];
V = g * cross(e3, q) - damp * w + k * (w'*q) * q;

F = zeros(6, 1);

F(1:3) = w;

F(4:6) = skw(q) * V;

end