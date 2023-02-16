clearvars
clc

Yn = normalize(rand(3,1));
my_zero = skw(zeros(3,1));
par = {[2 0 0; 0 2 0; 0 0 1];0};

my_x = expm(my_zero)*Yn;
my_field = field(my_x, par{1}, par{2});

rslt = dexpinv(my_zero, my_field);

disp(rslt)
disp(par{1}*(Yn*Yn')-Yn*Yn'*par{1})

function y = normalize(y)
    y = y/norm(y);
end

function A = field(x, D, C)
    f = (eye(3)-x*x')*D*x;
    A = x*f'-f*x'+C*skw(x);
end

function F = dexpinv(a,b)
    Theta = norm(a);
    theta = Theta * 0.5;
    aa = skw(a);
    if Theta > 1e-8
        fun = (1-theta*(cos(theta)/sin(theta)))/(Theta^2);
    else
        fun =  1/12 + Theta^2/720 + Theta^4/30240;
    end
    operator = eye(3) - 0.5*aa + fun*aa*aa;
    F = operator*b;
end