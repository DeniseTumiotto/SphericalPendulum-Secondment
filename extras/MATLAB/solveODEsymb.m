clearvars
close all
clc

syms u(t)

ode = diff(u) == -u-u.^3;

uSol(t) = dsolve(ode);
