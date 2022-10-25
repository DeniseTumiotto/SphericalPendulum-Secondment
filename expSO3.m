function rslt = expSO3(input)
%
% voglio trovare qual e` il vero spazio su cui lavoro
% e poi scrivere qui qual e` l'esponenziale del gruppo
%

[m, n] = size(input);

if m > 3 || n > 3
    w = input(1:3, end);
else
    w = input;
end

theta = norm(w);

rslt = eye(3) + sin(theta)/theta * skw(w) + (1-cos(theta))/theta^2 * skw(w)^2;

end