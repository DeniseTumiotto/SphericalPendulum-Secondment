function sol = LieEuler(vecField, action, exp, p, h)

% Computes one time step update with Lie Euler 

k1 = vecField(p);

sol = action(exp(h*k1), p); 

end