clearvars
close all
clc

syms A w q d
syms w1 w2 w3 q1 q2 q3

g = 9.81;
e3 = [0; 0; 1];

% w = [w1; w2; w3];
w = zeros(3,1);
q = [q1; q2; q3];

skw = @(x) [ 0     -x(3) x(2);
            x(3)  0     -x(1);
             -x(2) x(1)  0];

A = [skw(w) zeros(3); g.*skw(e3) -d*eye(3)];
% A = [zeros(3) eye(3); (norm(w)^2-q(3)*g)*eye(3) -d*eye(3)];
aux = [q*transpose(q)+(w*transpose(w))/(1+transpose(w)*w) (w*transpose(q))/(1+transpose(w)*w);
    (q*transpose(w))/(1+transpose(w)*w) (q*transpose(q))/(1+transpose(w)*w)];
J = eye(6)-aux;
A = J*A;

[T, lambda] = eig(A);

lambda = diag(lambda);
disp(lambda)

% C = inv(transpose(T))*inv(T);
% 
% assume(d>0)
% assume(w1>0)
% assume(w2>0)
% assume(w3>0)
% assume(q1>0)
% assume(q2>0)
% assume(q3<0)
% disp(sign(lambda))

% colorVec = ['k', 'k', 'k', 'k', 'b', 'r'];
% x = linspace(-1, 1, 10);
% 
% for n = 1:size(lambda, 1)
%     lambdaI = subs(lambda(n), {w1, w2, w3}, {x, x, x});
%     plot(real(lambdaI), imag(lambdaI), 'o', ...
%         'MarkerFaceColor', colorVec(n), 'MarkerSize', 5, ...
%         'MarkerEdgeColor', 'none');
%     hold on
% end
% 
% grid on
% hold off