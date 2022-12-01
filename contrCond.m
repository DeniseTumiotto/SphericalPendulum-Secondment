function [rslt, lam] = contrCond(space, y0, param)

f = vecField(space, y0, param);
J = jac(space, y0, param);

symmJ = symm_jac(f,J,y0);

lam = eig(symmJ);

if all(lam<0)
    rslt = 1;
else
    rslt = 0;
end

function this = symm_jac(A,B,y)
    this = [B(1,1),                          0.5*(B(1,2)/(cos(y(1))^2)+B(2,1));
            0.5*(cos(y(1))^2*B(2,1)+B(1,2)), B(2,2)-A(1)*tan(y(1))];
end
end

