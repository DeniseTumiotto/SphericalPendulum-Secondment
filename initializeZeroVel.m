function [q0, w0, z0, v] = initializeZeroVel()

% This method randomly picks a position point in S^2.
% The initial velocity remains zzero.

w0 = zeros(3, 1);
v = w0;
q0 = rand(3, 1);
q0(3) = -q0(3);
q0 = q0/norm(q0);

z0 = [q0; w0];

end