function rslt = actionSO3xR3(q, p)
%
% action on SO(3)xR3
%

v1 = q(1:3);
R1 = reshape(q(3 + (1:9)), 3, 3);
v2 = p(1:3);
R2 = reshape(p(3 + (1:9)), 3, 3);

rslt = zeros(size(q));

rslt(1:3) = v1 + v2;
rslt(3+(1:9)) = reshape(R1*R2, 9, 1);

end