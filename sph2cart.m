function rslt = sph2cart(q)
% Funtion from spherical coordinates to cartesian coordinate on S2

if size(q, 1) == 2 || size(q, 2) == 2
    q = [1, q(1), q(2)];
end

rslt = zeros(3,1);
rslt(1) = q(1)*cos(q(2))*cos(q(3));
rslt(2) = q(1)*cos(q(2))*sin(q(3));
rslt(3) = q(1)*sin(q(2));

end