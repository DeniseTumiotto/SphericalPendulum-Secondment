close all
clear
clc

data = readNPY('extras\test.npy');
[n1, n2, n3, n4] = size(data);

x1=reshape(data(1,1,:,:),n3,n4);
x2=reshape(data(1,2,:,:),n3,n4);
y1=reshape(data(2,1,:,:),n3,n4);
y2=reshape(data(2,2,:,:),n3,n4);
xx = linspace(0,2*pi,100);
x = cos(xx);
y = sin(xx);
z = zeros(1,100);

plotGradFlow(x1,x2)
plotGradFlow(y1,y2)
hold on
plot3(x,y,z)

rslt_expl = zeros(1,n4);
rslt_impl = zeros(1,n4);
rslt_expl_S = zeros(1,n4);
rslt_impl_S = zeros(1,n4);
rslt_expl_E = zeros(1,n4);
rslt_impl_E = zeros(1,n4);

nSteps = 5;
maxIt = 1000;
tol = 1e-6;
delta = 0.001;

for time = 1:n4
%     [v, w] = logMapS(x1(1:3, time), x1(4:6, time), ...
%                      x2(1:3, time), x2(4:6, time), ...
%                      nSteps, maxIt, tol, delta);
%     rslt_expl_S(time) = norm([v; w]);
% 
%     [v, w] = logMapS(y1(1:3, time), y1(4:6, time), ...
%                      y2(1:3, time), y2(4:6, time), ...
%                      nSteps, maxIt, tol, delta);
%     rslt_impl_S(time) = norm([v; w]);

    v = logMap(x1(1:3, time), x2(1:3, time));
    w = parallelTranslationAtoB(x2(1:3, time), x1(1:3, time), x2(4:6, time)) - x1(4:6, time);
    rslt_expl(time) = norm([v; w]);

    v = logMap(y1(1:3, time), y2(1:3, time));
    w = parallelTranslationAtoB(y2(1:3, time), y1(1:3, time), y2(4:6, time)) - y1(4:6, time);
    rslt_impl(time) = norm([v; w]);

%     rslt_expl_E(time) = euclidDistance(x1(:,time),x2(:,time));
%     rslt_impl_E(time) = euclidDistance(y1(:,time),y2(:,time));
end

figure()
plot(rslt_expl)
hold on
plot(rslt_impl)
% plot(rslt_expl_S)
% plot(rslt_impl_S)
% plot(rslt_expl_E)
% plot(rslt_impl_E)
legend('explicit distance', 'implicit distance')%, ...
%     'explicit distance Sasaki', 'implicit distance Sasaki', ...
%     'explicit distance Euclid', 'implicit distance Euclid')
grid on

figure()
plot(1:n4,x1(4,:), 1:n4,x1(5,:), 1:n4,x1(6,:))
grid on

figure()
plot(1:n4,x2(4,:), 1:n4,x2(5,:), 1:n4,x2(6,:))
grid on

figure()
plot(1:n4,y1(4,:), 1:n4,y1(5,:), 1:n4,y1(6,:))
grid on

figure()
plot(1:n4,y2(4,:), 1:n4,y2(5,:), 1:n4,y2(6,:))
grid on
