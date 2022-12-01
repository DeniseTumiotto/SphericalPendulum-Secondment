function f = vecField(space, y0, param)

if strcmp(space,'S2')
    cartcoo = sph2cart(y0);
    f = zeros(2,1);
    k = param.D(1,1)*cartcoo(1)^2+param.D(2,2)*cartcoo(2)^2+param.D(3,3)*cartcoo(3)^2;
    k1 = k-param.D(1,1);
    k2 = k-param.D(2,2);
    k3 = k-param.D(3,3);
    f(1) = k3*cartcoo(3)/sqrt(1-cartcoo(3)^2);
    f(2) = (k2-k1)*cartcoo(1)*cartcoo(2)/(cartcoo(1)^2+cartcoo(2)^2);
elseif strcmp(space,'TS2')
    f = zeros(4,1);

else
    disp('Not a valid option!')
end

end