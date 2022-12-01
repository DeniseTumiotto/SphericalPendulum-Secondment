function rslt = cart2sph(q)
% Funtion from cartesian coordinate on S2 to spherical coordinates

rslt(1) = sqrt(q(1)^2 + q(2)^2 + q(3)^2);
rslt(2) = atan(q(3)/sqrt(q(1)^2+q(2)^2));
rslt(3) = mod(atan(q(2)/q(1)),2*pi);

end