function rslt = forwardEuler(y, f, dt)

rslt = y + dt * f(y);

end
