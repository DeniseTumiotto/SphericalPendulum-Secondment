function sol = NewtonIt(residual, jacobian, p, h, max_it, atol, rtol)

% Newton iteration for implicit methods
% residual is a function to evaluate the RHS of the system
% jacobian evaluates the iteration matrix to perform the Newton iteration
% p is the initial value of the current iteration
% h is the time step size
% max_it is the maximum number of iteration
% tol is the tolerance to determine when to stop the integration step

pNew = p;
dx = zeros(size(p));
for i = 1:max_it
    sol = residual(p, pNew, h);

    res = -pNew + sol;
    if norm(res) < rtol*norm(dx) + atol
%             disp(i)
        break
    end

    S_t = jacobian(pNew, p, h);
    dx = -S_t\res;

    pNew = pNew + dx;
end

if i == max_it
    disp('not converged!')
end

end