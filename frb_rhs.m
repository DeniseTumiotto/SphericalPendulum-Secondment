function rslt = frb_rhs(yOld, yNew, h, field, action, exponential, mode)
%
% Descrete right hand side for implicit method
% function for the residual
%

rslt = zeros(6, 1);

switch mode
    case "implicit Lie Euler method"
        rslt = action(exponential(h*field(yNew)), yOld);
    case "implicit midpoint rule"
        rslt = action(exponential(h*field((yNew+yOld)/2)), yOld);

end