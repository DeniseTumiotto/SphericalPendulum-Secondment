function J = jac(space, y0, param)

if strcmp(space,'S2')
    % ATTENZIONE! sto usando il fatto che D(1,1)=D(2,2)
    J = [(param.D(1,1)-param.D(3,3)) * ( cos(y0(1))^2 - sin(y0(1))^2 * cos(y0(1))/(abs(cos(y0(1)))) ) 0;
        0 (param.D(1,1)-param.D(3,3)) * ( - sin(y0(1))^2 * cos(y0(1))/(abs(cos(y0(1)))) )];
elseif strcmp(space,'TS2')
    J = zeros(4);
else
    disp('Not a valid option!')
end

end