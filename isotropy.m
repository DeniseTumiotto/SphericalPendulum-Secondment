clc

full = 0;

% reading solution with different values of C
s = readpy('isotropy',1,1);
n = size(s{1});

% copying (manually) c vector
C = [0, 2, 5, 10, 20];

color = linspecer(n(1)+2);

% reproduce graph on the sphere of the solution

figure()
hold on
% plot data
for k = 1:n(1)
    x = reshape(s{1}(k,:,:),n(2),n(3));
    radius = sqrt(x(1,1)^2+x(2,1)^2);
    % plot3(radius * cos(linspace(0,2*pi,100)), radius * sin(linspace(0,2*pi,100)), x(3,1) * ones(100,1),'-k','LineWidth',1)
    plot3(x(1,:),x(2,:),x(3,:), "-", 'LineWidth', 3,'Color',color(k,:),'DisplayName',strcat('c=',num2str(C(k))))
end
legend('FontSize',18,'AutoUpdate','off')

% Create sphere surface
[xS2, yS2, zS2] = sphere(360);
small = 360/2;
if ~full
    % only lower half of the sphere
    h = surf(xS2(1:small,:), yS2(1:small,:), zS2(1:small,:), 'FaceAlpha', 0.5); 
    h.EdgeColor = 'none';
    plot3(xS2(small,:),yS2(small,:),zS2(small,:),'-k',LineWidth=1)
    % colormap('bone')
    % colormap('pink')
else
    % complete sphere
    h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
    h.EdgeColor = 'none';
end

axis equal


