function plotGradFlow(y1,y2)
% plots the solutions of a gradient flow on S2
% 
% :param y1: solution 1
% :param y2: solution 2
%
% :returns: plot the two solution on the unitary sphere
%


s0(1,:) = cart2sph(y1(:,1));
s0(2,:) = cart2sph(y2(:,1));

figure()

% Create sphere surface
[xS2, yS2, zS2] = sphere(360);
h = surf(xS2(1:180,:), yS2(1:180,:), zS2(1:180,:), 'FaceAlpha', 0.1); 
h.EdgeColor = 'none';
hold on
plot3(cos(0:0.001:2*pi),sin(0:0.001:2*pi),zeros(size(0:0.001:2*pi)),'-k',LineWidth=1)
% Plot initial values
plot3(y1(1,1), y1(2,1), y1(3,1), 'o', 'MarkerSize', 5, ...
'MarkerFaceColor', 'green', ...
'MarkerEdgeColor', 'none')
plot3(y2(1,1), y2(2,1), y2(3,1), 'o', 'MarkerSize', 5, ...
'MarkerFaceColor', 'red', ...
'MarkerEdgeColor', 'none')
plot3(y1(1,:),y1(2,:),y1(3,:), "-o", 'MarkerSize', 3, ...
'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', ...
'LineWidth', 3,'Color','r')
plot3(y2(1,:),y2(2,:),y2(3,:), "-o", 'MarkerSize', 3, ...
'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', ...
'LineWidth', 3,'Color','b')


title({['($\phi_0,\theta_0$)=(',num2str(s0(1,2)),', ',num2str(s0(1,3)),'), '], ...
       ['($\tilde{\phi}_0,\tilde{\theta}_0$)=(',num2str(s0(2,2)),', ',num2str(s0(2,3)),')']}, ...
        'Interpreter', 'latex','FontSize',24)
xlabel('x',FontSize=18)
ylabel('y',FontSize=18)
zlabel('z',FontSize=18)
end