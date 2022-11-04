close all
clear
clc

% exact sol contractive, implicit no
% s0(1,:) = [-0.19194617698007518,5.988667862558153];
% s0(2,:) = [-0.524058693543719,2.0227585780325947];
s0(1,:) = [0.48149073285181876,3.225938294516135];
s0(2,:) = [1.1559010717705054,2.2361071703320317];

% exact sol NOT contractive, implicit yes
% s0(1,:) = [1.3780657846808189,3.796433641270566];
% s0(2,:) = [-0.5247740271925787,0.1765692724968637];
% s0(1,:) = [-0.0056398875912720015,5.0894692582092125];
% s0(2,:) = [-0.19000492393810275,1.2795274856218333];
% s0(1,:) = [1.0686561239381236,0.1006863545605521];
% s0(2,:) = [0.9105258034060973,5.366661364741157];

% both contractive
% s0(1,:) = [-0.5118121864624298,5.593179094130903];
% s0(2,:) = [0.06437717050937097,3.3171602172857084];

q(:,1) = sph2cart(s0(1, :));
q(:,2) = sph2cart(s0(2, :));

nSpace = 1000;
t = linspace(0,1,nSpace);
geo = zeros(3,nSpace);
k = 0;
for i = t
    k = k + 1;
    geo(:,k) = ((1-i).*q(:,1)+i.*q(:,2))/norm((1-i).*q(:,1)+i.*q(:,2));
end

figure()
myColor = colormap("hsv");
[m, ~] = size(myColor);

% Create sphere surface
[xS2, yS2, zS2] = sphere(360);
h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
h.EdgeColor = 'none';
hold on
plot3(q(1, 1), q(2, 1), q(3, 1), 'o', 'MarkerSize', 5, ...
'MarkerFaceColor', 'green', ...
'MarkerEdgeColor', 'none')
plot3(q(1, 2), q(2, 2), q(3, 2), 'o', 'MarkerSize', 5, ...
'MarkerFaceColor', 'blue', ...
'MarkerEdgeColor', 'none')

plot3(geo(1,:), geo(2,:), geo(3,:), 'LineWidth', 5, 'Color', 'r')
plot3(1,0,0, 'ok', 'MarkerSize',5, 'MarkerFaceColor', 'k')


qe1 = zeros(3,nSpace);
qe2 = zeros(3,nSpace);
k = 0;
for i = linspace(0.01,15,nSpace)
    k = k+1;
    qe1(:,k) = exactSol(s0(1,:),i);
    qe2(:,k) = exactSol(s0(2,:),i);
    plot3(qe1(1,k),qe1(2,k),qe1(3,k), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', 'green', ...
    'MarkerEdgeColor', 'none')
    plot3(qe2(1,k),qe2(2,k),qe2(3,k), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', 'blue', ...
    'MarkerEdgeColor', 'none')
end

% plot3(0.48882947,0.04104581,0.87141321, 'xk', 'MarkerSize', 5, 'LineWidth',12)
% plot3(0.44838799,-0.42876861,0.78428674, 'xk', 'MarkerSize', 5, 'LineWidth',12)
% plot3(-0.73603298, -0.4975225, 0.45904991, 'xk', 'MarkerSize', 5, 'LineWidth',12)
% plot3(-0.11765011, 0.39620143, 0.91059479, 'xk', 'MarkerSize', 5, 'LineWidth',12)
% plot3(0.73378779, -0.47539177, -0.48534333, 'xk', 'MarkerSize', 5, 'LineWidth',12)
% plot3(-0.78228239, -0.61965617, 0.06372197, 'xk', 'MarkerSize', 5, 'LineWidth',12)

title(['($\phi_0,\theta_0$)=(',num2str(s0(1,1)),', ',num2str(s0(1,2)),'), ', ...
       '($\tilde{\phi}_0,\tilde{\theta}_0$)=(',num2str(s0(2,1)),', ',num2str(s0(2,2)),')'], ...
        Interpreter='latex')

function rslt = exactSol(x0, t)
    phi = exp(-0.1*t)*x0(1);
    theta = exp(-2*t)*x0(2);
    
    rslt = sph2cart([phi,theta]);
end
