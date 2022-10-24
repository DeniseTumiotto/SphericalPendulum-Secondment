function animationSphere(sols, params, step)

if nargin < 3
    step = 1;
end

[n, ~] = size(sols);
color = ['r', 'b'];

figure('Units','centimeters','Position',[10 10 20 15])
% Create sphere surface
[xS2, yS2, zS2] = sphere(360);
h = surf(xS2, yS2, zS2, 'FaceAlpha', 0.1); 
h.EdgeColor = 'none';
hold on
for m = 1:n
    q = sols{m}(1:3, :);
    steps = params{m}.N_TIME;
    dstep = 1 / steps;
    for i = 1:step:steps
        % Plot pendulum at time step i
        plot3([0, q(1,i)], [0, q(2,i)], [0, q(3,i)], '*', 'Color', color(m));
        xlabel("x")
        ylabel("y")
        zlabel("z")
        str = "Time evolution of the pendulums, T="+num2str(dstep*i);
        axis equal;
        title(str, 'FontSize', 16)
        pause(0.0000000000001);
    end
end
end