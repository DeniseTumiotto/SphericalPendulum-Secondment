syms x y z
syms q1 q2 q3
syms a c

A = [0 -z y;
    z 0 -x;
    -y x 0];

q = [q1; q2; q3];

b = [(c-a) * q(3)^2 * q(1);
    (c-a) * q(3)^2 * q(2);
    (c-a) * (q(3)^2-1) * q(3)];


eqns = A*q-b==0;

S = solve(eqns,[x y z]);