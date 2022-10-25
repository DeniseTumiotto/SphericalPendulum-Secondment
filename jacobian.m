function sol = jacobian(v0, v, h, f, action, exponential, method)

% Jacobian of the residual value
% evaluation with Finite Differences

[m, ~] = size(v);

sol = zeros(m);
sol0 = action(exponential(h*f(v0)), v0);

switch method
    
    % Implicit Lie Euler
    case "implicit Lie Euler method"
        for i = 1:m
            dx = v(i)-v0(i);
            if dx ~= 0
                newV = v0;
                newV(i) = v(i);
                k1 = f(newV);
                partialSol = action(exponential(h*k1), v0);
                sol(:,i) = (partialSol-sol0)/dx;
            end
        end
        
        sol = sol-eye(m);

    % Implicit Midpoint Rule
    case "implicit midpoint rule"
        for i = 1:m
            dx = v(i)-v0(i);
            if dx ~= 0
                newV = v0;
                newV(i) = v(i);
                k1 = f((newV+v0)/2);
                partialSol = action(exponential(h*k1), v0);
                sol(:,i) = (partialSol-sol0)/(2*dx);
            end
        end
        
        sol = sol-eye(m);

    otherwise
        error('Method not implemented!')

end

end