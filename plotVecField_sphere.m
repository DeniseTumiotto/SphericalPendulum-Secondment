close all
% Energy and Gradient Functions
E = @(q,D) 0.5*q.'*D*q;
grad = @(q,D,f) - D*q + 2*f(q,D)*q;

% D MATRIX
a = 0.5;
c = 0.25;
inertia = [a 0 0; 0 a 0; 0 0 c];

% Color
r = 255/255;
g = 253/255;
b = 208/255;

% Create sphere surface
figure()
[xS2, yS2, zS2] = sphere(359);
h = surf(xS2, yS2, zS2); 
h.EdgeColor = 'none';
h.FaceColor = [r,g,b];
hold on
% plot3(cos(0:0.001:2*pi),sin(0:0.001:2*pi),zeros(size(0:0.001:2*pi)),'-k',LineWidth=1)
% plot3(0,0,1,'xk','MarkerSize',5)
% plot3(0,0,-1,'xk','MarkerSize',5)

N = 29;
[xS2, yS2, zS2] = sphere(N);
V = zeros(N+1,N+1);
U = zeros(N+1,N+1);
W = zeros(N+1,N+1);
for t = 1:N+1
    for j = 1:N+1
        tmpvel = grad([xS2(t,j);yS2(t,j);zS2(t,j)],inertia,E);
        V(t,j) = tmpvel(1);
        U(t,j) = tmpvel(2);
        W(t,j) = tmpvel(3);
    end
end

% Color
r = 51/255;
g = 51/255;
b = 178/255;

quiver3(xS2(1:floor((N+1)/2),:),yS2(1:floor((N+1)/2),:),zS2(1:floor((N+1)/2),:),...
    V(1:floor((N+1)/2),:),U(1:floor((N+1)/2),:),W(1:floor((N+1)/2),:),'Color',[r,g,b])
quiver3(xS2(ceil((N+1)/2)+1:end,:),yS2(ceil((N+1)/2)+1:end,:),zS2(ceil((N+1)/2)+1:end,:),...
    V(ceil((N+1)/2)+1:end,:),U(ceil((N+1)/2)+1:end,:),W(ceil((N+1)/2)+1:end,:),'Color',[r,g,b])
axis equal
set(gca,'XColor','none','YColor','none','ZColor','none')
set(gca, 'color', 'none')
grid('off')

% export_fig test.png -m3 -transparent