function plotGradFlow(y, small, separate, full)
% plots the solutions of a gradient flow on S2
% 
% :param y1: solution 1
% :param y2: solution 2
%
% :returns: plot the two solution on the unitary sphere
%
if nargin < 4
    full = 0;
    if nargin < 3
        separate = 0;
        if nargin < 2
            small = 360/2;
        end
    end
end

n = size(y,1);
s0 = cell(n,1);
color = linspecer(n+2);

if ~separate
    figure()
    % Create sphere surface
    [xS2, yS2, zS2] = sphere(360);
    if ~full
        h = surf(xS2(1:small,:), yS2(1:small,:), zS2(1:small,:), 'FaceAlpha', 0.1); 
        h.EdgeColor = 'none';
        hold on
        plot3(xS2(small,:),yS2(small,:),zS2(small,:),'-k',LineWidth=1)
    else
        h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
        h.EdgeColor = 'none';
        hold on
    end
end

for i = 1:n
    s0{i} = cart2sph(y{i}(:,1));
    
    if separate
        figure()
        % Create sphere surface
        [xS2, yS2, zS2] = sphere(360);
        if ~full
            h = surf(xS2(1:small,:), yS2(1:small,:), zS2(1:small,:), 'FaceAlpha', 0.1); 
            h.EdgeColor = 'none';
            hold on
            plot3(xS2(small,:),yS2(small,:),zS2(small,:),'-k',LineWidth=1)
        else
            h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
            h.EdgeColor = 'none';
            hold on
        end
    end
    
    % Plot initial values
    plot3(y{i}(1,1), y{i}(2,1), y{i}(3,1), 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'green', ...
    'MarkerEdgeColor', 'none')
    radius = sqrt(y{i}(1,1)^2+y{i}(2,1)^2);
    plot3(radius * cos(linspace(0,2*pi,100)), radius * sin(linspace(0,2*pi,100)), y{i}(3,1) * ones(100,1),'-k','LineWidth',1)
    plot3(y{i}(1,:),y{i}(2,:),y{i}(3,:), "-o", 'MarkerSize', 3, ...
    'MarkerEdgeColor', 'none', ...
    'LineWidth', 3,'Color',color(i,:))    
end
plot3(0,0,-1,'ok','MarkerFaceColor','k','MarkerSize',4)
plot3(0,0,1,'ok','MarkerFaceColor','k','MarkerSize',4)
title('Trajectories for c=0','FontSize',24)
xlabel('x',FontSize=20)
ylabel('y',FontSize=20)
zlabel('z',FontSize=20)
axis equal
end