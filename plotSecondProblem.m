function plotSecondProblem(y1,y2)

[~, nTime] = size(y1);

s0(1,:) = cart2sph(y1(:,1));
s0(2,:) = cart2sph(y2(:,1));

nSpace = 1000;
t = linspace(0,1,nSpace);
geo = zeros(3,nSpace);
k = 0;
for i = t
    k = k + 1;
    geo(:,k) = ((1-i).*y1(:,1)+i.*y2(:,1))/norm((1-i).*y1(:,1)+i.*y2(:,1));
end

figure()

% Create sphere surface
[xS2, yS2, zS2] = sphere(360);
h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
h.EdgeColor = 'none';
hold on
% Plot initial values
plot3(y1(1,1), y1(2,1), y1(3,1), 'o', 'MarkerSize', 5, ...
'MarkerFaceColor', 'green', ...
'MarkerEdgeColor', 'none')
plot3(y2(1,1), y2(2,1), y2(3,1), 'o', 'MarkerSize', 5, ...
'MarkerFaceColor', 'red', ...
'MarkerEdgeColor', 'none')
% Plot geodesic between initial values and equilibrium point
plot3(geo(1,:), geo(2,:), geo(3,:), 'LineWidth', 5, 'Color', 'b')
plot3(1,0,0, 'ok', 'MarkerSize',5, 'MarkerFaceColor', 'k')

% Evaluate and plot exact solution
qe1 = zeros(3,nSpace);
qe2 = zeros(3,nSpace);
k = 0;
for i = linspace(0.01,15,nSpace)
    k = k+1;
    qe1(:,k) = exactSol(s0(1,2:3),i);
    qe2(:,k) = exactSol(s0(2,2:3),i);
    plot3(qe1(1,k),qe1(2,k),qe1(3,k), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', 'green', ...
    'MarkerEdgeColor', 'none')
    plot3(qe2(1,k),qe2(2,k),qe2(3,k), 'o', 'MarkerSize', 3, ...
    'MarkerFaceColor', 'red', ...
    'MarkerEdgeColor', 'none')
end

for k = 1:nTime
    plot3(y1(1,k),y1(2,k),y1(3,k), "pentagram", 'MarkerSize', 3, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'none')
    plot3(y2(1,k),y2(2,k),y2(3,k), "pentagram", 'MarkerSize', 3, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'none')
end

title(['($\phi_0,\theta_0$)=(',num2str(s0(1,2)),', ',num2str(s0(1,3)),'), ', ...
       '($\tilde{\phi}_0,\tilde{\theta}_0$)=(',num2str(s0(2,2)),', ',num2str(s0(2,3)),')'], ...
        'Interpreter', 'latex')

function rslt = exactSol(x0, t)
    phi = exp(-0.1*t)*x0(1);
    theta = exp(-2*t)*x0(2);
    
    rslt = sph2cart([phi,theta]);
end
end