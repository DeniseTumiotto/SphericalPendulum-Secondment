function plotGradFlow(y1,y2)

[~, nTime] = size(y1);

s0(1,:) = cart2sph(y1(:,1));
s0(2,:) = cart2sph(y2(:,1));

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

for k = 1:nTime
    plot3(y1(1,k),y1(2,k),y1(3,k), "o", 'MarkerSize', 3, ...
    'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'none')
    plot3(y2(1,k),y2(2,k),y2(3,k), "o", 'MarkerSize', 3, ...
    'MarkerFaceColor', 'b', ...
    'MarkerEdgeColor', 'none')
end

title(['($\phi_0,\theta_0$)=(',num2str(s0(1,2)),', ',num2str(s0(1,3)),'), ', ...
       '($\tilde{\phi}_0,\tilde{\theta}_0$)=(',num2str(s0(2,2)),', ',num2str(s0(2,3)),')'], ...
        'Interpreter', 'latex')
end