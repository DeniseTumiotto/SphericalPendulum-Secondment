function rslt = expSO3xR3(input)
%
% exponential map in SO3xR3
%

% rslt = zeros(3, 4);

A = expSO3(input);

rslt = [A; input(1:3)];


end